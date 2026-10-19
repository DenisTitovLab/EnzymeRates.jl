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
    analytical_peak_turnover_fn::Union{Function,Nothing} = nothing

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
    raw::Union{EnzymeRates.RawMechanism,EnzymeRates.RawAllostericMechanism}
end

"""
    _testhelper_mechanism_call(call, who) -> (allosteric, raw)

Whether `call`, an `@enzyme_mechanism` or `@allosteric_mechanism` call, is the allosteric
one, and the raw value the package's DSL walker reads from its block; an error naming the
test macro `who` for any other expression.
"""
function _testhelper_mechanism_call(call, who::String)
    head = call isa Expr && call.head === :macrocall ? call.args[1] : nothing
    head in (Symbol("@enzyme_mechanism"), Symbol("@allosteric_mechanism")) ||
        error("$who wraps an @enzyme_mechanism or @allosteric_mechanism call; got $call")
    allosteric = head === Symbol("@allosteric_mechanism")
    allosteric, EnzymeRates._raw_mechanism(call.args[end], allosteric)
end

"""
    @_testhelper_with_raw call

`(mechanism, raw)` for `call`, an `@enzyme_mechanism` or `@allosteric_mechanism` call:
the compiled mechanism the call builds, and the raw value the package's DSL walker
reads from its block, which keeps the steps and regulatory sites as written. It fills
`_testhelper_MechanismTestSpec.raw`.
"""
macro _testhelper_with_raw(call)
    _, raw = _testhelper_mechanism_call(call, "@_testhelper_with_raw")
    esc(:(($call, $raw)))
end

# Every declared substrate and product appears in some step: some step takes it up or
# gives it off. The constructor rejects a mechanism that declares a reactant no step
# exchanges, so this checks that the test helpers' reader (`_testhelper_graph`) reads
# every exchange.
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

"""Whether the concrete mechanism `x` reads back from the text `show` prints for it
(`_parse_mechanism`) as itself."""
_testhelper_reads_back(x) = EnzymeRates._parse_mechanism(sprint(show, x)) == x

"""
    _testhelper_graph(m) -> NamedTuple

`m` (concrete or compiled) as the test helpers read it, a `NamedTuple`: `forms`, each
form's name (`Vector{String}`, stored order); `free`, whether each form holds no
metabolite and no residual; `steps`, one `NamedTuple` per step, group by group, each group
in edge order, in stored direction, with `from` and `to` (form positions), `consumed`
and `released` (metabolite names, a copy by its reactant's name), `consumed_modifier`
and `released_modifier` (which of them bind as modifiers), `rapid_equilibrium`, `group`,
and `names` and `inactive_names`, the package names of its constants in the active (or
only) and the inactive conformation (empty for a step the inactive law does not keep);
`tags`; `inactive`, which steps the inactive law keeps;
`catalytic_multiplicity`; `sites`, each `(; multiplicity, ligands, allo_states)`;
`substrates`, `products` and `modifiers`. A plain mechanism has no tags, no inactive
conformation (`inactive` and every `inactive_names` empty), no sites and
`catalytic_multiplicity` `nothing`. Forms come in the canonical form order: free E first,
then by content.
"""
function _testhelper_graph(m)
    mech = _testhelper_lift(m)
    allo = mech isa EnzymeRates.AllostericMechanism
    cm = EnzymeRates.catalytic(mech)
    kept = EnzymeRates._kept_members(mech)
    tags = allo ? copy(EnzymeRates.cat_allo_states(mech)) : Symbol[]
    # The edges group by group, each group's in edge order, as the derivation reads them.
    edges = EnzymeRates._flat_edges(cm)
    _testhelper_names(k, state) = Symbol[EnzymeRates.name(p, mech, kept)
                                         for p in EnzymeRates._edge_constants(cm, k, state)]
    in_law = allo ? EnzymeRates._inactive_mask(mech) : falses(length(cm.edges))
    steps = map(edges) do k
        e, g = cm.edges[k], cm.group[k]
        consumed = cm.table[EnzymeRates.entries(e.consumed)]
        released = cm.table[EnzymeRates.entries(e.released)]
        shared = allo && tags[g] === :EqualAI
        (from = e.from, to = e.to,
         consumed = Symbol[x.name for x in consumed],
         released = Symbol[x.name for x in released],
         consumed_modifier = Bool[x.modifier for x in consumed],
         released_modifier = Bool[x.modifier for x in released],
         rapid_equilibrium = e.rapid_equilibrium, group = g,
         names = _testhelper_names(k, !allo ? :None : shared ? :EqualAI : :A),
         inactive_names = in_law[k] ? _testhelper_names(k, shared ? :EqualAI : :I) :
                          Symbol[])
    end
    sites = [(multiplicity = EnzymeRates.multiplicity(site),
              ligands = copy(EnzymeRates.ligands(site)),
              allo_states = copy(EnzymeRates.allo_states(site)))
             for site in (allo ? EnzymeRates.regulatory_sites(mech) :
                                 EnzymeRates.RegulatorySite[])]
    (forms = String[EnzymeRates._form_name(cm, f) for f in eachindex(cm.forms)],
     free = Bool[EnzymeRates._is_free(f) for f in cm.forms],
     steps, tags,
     inactive = Int[j for (j, k) in enumerate(edges) if in_law[k]],
     catalytic_multiplicity = allo ? EnzymeRates.catalytic_multiplicity(mech) : nothing,
     sites, substrates = EnzymeRates.substrates(mech),
     products = EnzymeRates.products(mech), modifiers = EnzymeRates.modifiers(mech))
end

"""The name of every step constant of the plain mechanism `m`, fitted or dependent: each
step's constants as the naming chokepoint renders them."""
_testhelper_step_constant_names(m) =
    Set(x for s in _testhelper_graph(m).steps for x in s.names)

# The names of a compiled mechanism's substrates and products, and of its regulators: the
# modifiers an EnzymeMechanism's steps bind.
_testhelper_substrates(m) = _testhelper_graph(m).substrates
_testhelper_products(m) = _testhelper_graph(m).products
_testhelper_regulators(m::EnzymeRates.EnzymeMechanism) = _testhelper_graph(m).modifiers

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
`m`, compiled or concrete, with one catalytic subunit: an allosteric mechanism rebuilt
with catalytic multiplicity 1, whose MWC law is the exact steady state of one protomer; a
plain mechanism unchanged.
"""
_testhelper_single_subunit(m::EnzymeRates.EnzymeMechanism) = m
_testhelper_single_subunit(m::EnzymeRates.AllostericEnzymeMechanism) =
    EnzymeRates.AllostericEnzymeMechanism(
        _testhelper_single_subunit(EnzymeRates.AllostericMechanism(m)))
_testhelper_single_subunit(m::EnzymeRates.Mechanism) = m
_testhelper_single_subunit(am::EnzymeRates.AllostericMechanism) =
    EnzymeRates.AllostericMechanism(am.catalytic, copy(am.cat_allo_states), 1,
                                    am.regulatory_sites)

"""Whether `s`, a step of `_testhelper_graph`, is a binding: it takes up one metabolite and
gives off none."""
_testhelper_is_binding(s) = length(s.consumed) == 1 && isempty(s.released)

"""Whether `s`, a step of `_testhelper_graph`, is an isomerization: it takes up and gives
off nothing."""
_testhelper_is_iso(s) = isempty(s.consumed) && isempty(s.released)

"""The steps of `g` (`_testhelper_graph`) kinetic group by kinetic group, in group
order."""
_testhelper_groups(g) = [[s for s in g.steps if s.group == k]
                         for k in unique(s.group for s in g.steps)]

"""The first step of each kinetic group of `g` (`_testhelper_graph`), in group order."""
_testhelper_first_steps(g) = first.(_testhelper_groups(g))

"""
    _testhelper_re_segments(g, kept = trues(length(g.steps))) -> (segments, extras)

The rapid-equilibrium segments of `g` (`_testhelper_graph`) over the steps `kept` keeps,
read off its steps alone: `segments`, each segment's forms breadth-first from its first
form (each form's rapid-equilibrium steps in step order), the segments in the order of
their first forms; `extras`, for each form some kept step touches, in form order, how many
more of each metabolite name it carries than the lowest form of its segment, the names it
carries no more of left out. Rapid equilibrium fixes the ratio of two forms of a segment
as a monomial, one factor per metabolite bound or released between them.
"""
function _testhelper_re_segments(g, kept = trues(length(g.steps)))
    nf = length(g.forms)
    present = falses(nf)
    adjacent = [Tuple{Int, Dict{Symbol, Int}}[] for _ in 1:nf]
    for (j, s) in enumerate(g.steps)
        kept[j] || continue
        present[s.from] = present[s.to] = true
        s.rapid_equilibrium || continue
        change = Dict{Symbol, Int}()
        foreach(x -> change[x] = get(change, x, 0) + 1, s.consumed)
        foreach(x -> change[x] = get(change, x, 0) - 1, s.released)
        push!(adjacent[s.from], (s.to, change))
        push!(adjacent[s.to], (s.from, Dict(x => -n for (x, n) in change)))
    end
    exponent = [Dict{Symbol, Int}() for _ in 1:nf]
    reached = falses(nf)
    segments = Vector{Int}[]
    for root in findall(present)
        reached[root] && continue
        reached[root] = true
        segment = [root]
        i = 1
        while i <= length(segment)
            u = segment[i]
            i += 1
            for (v, change) in adjacent[u]
                reached[v] && continue
                reached[v] = true
                exponent[v] = mergewith(+, exponent[u], change)
                push!(segment, v)
            end
        end
        for x in union(Set{Symbol}(), (keys(exponent[f]) for f in segment)...)
            lowest = minimum(get(exponent[f], x, 0) for f in segment)
            foreach(f -> exponent[f][x] = get(exponent[f], x, 0) - lowest, segment)
        end
        push!(segments, segment)
    end
    segments, [filter(p -> p.second > 0, exponent[f]) for f in findall(present)]
end

"""The text of `s`, a step of `g` (`_testhelper_graph`), each side its form's name and its
metabolites' names, a modifier copy by its reactant's name: `"E + S ⇌ ES"`,
`"ES <--> EP"`, `"ES + P ⇌ EPmodS"`."""
_testhelper_step_text(g, s) =
    join([g.forms[s.from]; String.(s.consumed)], " + ") *
    (s.rapid_equilibrium ? " ⇌ " : " <--> ") *
    join([g.forms[s.to]; String.(s.released)], " + ")

"""Edge `k` of the catalytic graph of `m` (concrete or compiled) as text, in its stored
direction, each side its form's name (`_form_name`) and the labels of the metabolites the
edge takes up or gives off there (`_entry_label`, so a reactant's modifier copy reads
`Amod`): `"E + S ⇌ ES"`, `"EA + B <--> EQ + P"`, `"E + Amod ⇌ EAmod"`."""
function _testhelper_edge_text(m, k::Int)
    g = EnzymeRates.catalytic(_testhelper_lift(m))
    e = g.edges[k]
    labels(c) = [EnzymeRates._entry_label(g, i) for i in EnzymeRates.entries(c)]
    side(f, c) = join([EnzymeRates._form_name(g, f); labels(c)], " + ")
    side(e.from, e.consumed) * (e.rapid_equilibrium ? " ⇌ " : " <--> ") *
        side(e.to, e.released)
end

"""The edge of the catalytic graph of `m` whose text (`_testhelper_edge_text`) is `text`."""
_testhelper_edge(m, text) =
    only(k for k in eachindex(EnzymeRates.catalytic(_testhelper_lift(m)).edges)
         if _testhelper_edge_text(m, k) == text)

"""
    _testhelper_consistent_point(m, rng) -> (; vals, Keq, active, inactive, L, K, kf)

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

`K[k, c]` is edge `k`'s equilibrium constant in its stored direction in conformation `c`
and `kf[k, c]`, for a steady-state edge, its forward rate, for every edge of the
catalytic graph, the ones the inactive law deletes or prunes included; a plain
mechanism repeats column 1 in column 2. The law's constants are among them.
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
    # Every edge's equilibrium constant in both conformations and, for a steady-state
    # edge, its forward rate: the law's own where the law holds the edge, an `:EqualAI`
    # edge's active one, and otherwise one draw per kinetic group and conformation (each
    # member of an `:OnlyA` group in the inactive conformation its own). The
    # well-posedness oracles read these for the edges the inactive law deletes or prunes.
    # (`forward` holds the forward rates: a local named `kf` would be the one the closure
    # `_testhelper_conformation_steps` assigns.)
    cm = EnzymeRates.catalytic(mech)
    flat = EnzymeRates._flat_edges(cm)
    K = ones(BigFloat, length(cm.edges), 2)
    forward = ones(BigFloat, length(cm.edges), 2)
    rates = Dict{Tuple{Int, Int}, BigFloat}()
    for c in 1:nconf, s in g.steps
        names = c == 1 || tag(s.group) === :EqualAI ? s.names : s.inactive_names
        s.rapid_equilibrium || isempty(names) || (rates[(s.group, c)] = vals[names[1]])
    end
    for c in 1:nconf, (j, s) in enumerate(g.steps)
        k = flat[j]
        K[k, c] = exp(u[(c - 1) * nf + s.from] + sum(gm, s.consumed; init = big(0)) -
                      u[(c - 1) * nf + s.to] - sum(gm, s.released; init = big(0)))
        s.rapid_equilibrium && continue
        key = c == 2 && tag(s.group) === :OnlyA ? (-k, 2) : (s.group, c)
        forward[k, c] = get!(() -> exp(4 * big(rand(rng)) - 2), rates, key)
    end
    allo || (K[:, 2] .= K[:, 1]; forward[:, 2] .= forward[:, 1])
    (; vals, Keq, active, inactive, L, K, kf = forward)
end

"""
    _testhelper_steady_state_rate(m, point, concs; Λ = big(10)^30, flip = 1) -> BigFloat

The exact steady-state rate per enzyme of the mass-action network of `m` at the constants
of `point` (`_testhelper_consistent_point`) and the concentrations `concs`, in 256-bit
BigFloat: the net consumption of the first substrate. Each rapid-equilibrium step runs
as a fast pair at rate scale `Λ`, so the steady state approaches the rapid-equilibrium
limit with a relative error of about 1/Λ times the spread of the constants. An allosteric
mechanism (one catalytic subunit) adds the inactive conformation's forms and steps and one
flip between the two copies of form `flip`, free E unless a test asks for another, at
equilibrium ratio `L` times the ratio of the conformations' regulatory-site factors. Errors
when the network is disconnected, since its steady state is then not unique.
"""
function _testhelper_steady_state_rate(m, point, concs::AbstractDict; Λ = big(10)^30,
                                       flip::Int = 1)
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
            push!(edges, (id(:A, flip), id(:I, flip),
                          _testhelper_flip_ratio(mech, point, concs), big(1), 0))
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
The rate law of `m` per enzyme, evaluated in 256-bit BigFloat from the polynomials
`rate_equation` compiles (`_num_den_exprs`) at the parameter values `params` (by name: the
fitted constants, `Keq` and, for an allosteric mechanism, `L` and the regulator
constants) and the concentrations `concs`. A compiled `m` goes through its derivation
entry, which runs the boundary check; a concrete `m` is read as it is, so the oracles can
evaluate the law of a mechanism the guard rejects, though class 5 still refuses a law that
is 0/0 on a face (`_num_den_exprs` throws its message). The dependent constants are
computed from `_dependent_param_exprs`, as the compiled body computes them.
"""
function _testhelper_law_value(m, params::AbstractDict, concs::AbstractDict)
    setprecision(BigFloat, 256) do
        mech = m isa EnzymeRates.AbstractEnzymeMechanism ?
               EnzymeRates._derivation_input(typeof(m)) : m
        num, den = EnzymeRates._num_den_exprs(mech)
        dep, _ = EnzymeRates._dependent_param_exprs(mech)
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
    _testhelper_law_inputs(m, fitted, rng) -> (point, concs, params)

What a check of a rate law against an oracle reads, drawn from `rng` in this order: a
consistent point of the concrete mechanism `m` at 256 bits (`_testhelper_consistent_point`),
then each metabolite's concentration, log-uniform between e^-1.5 and e^1.5, in
`metabolites(m)` order. `params` holds the point's values of the fitted constants `fitted`
and its `Keq`, the law's inputs.
"""
function _testhelper_law_inputs(m, fitted, rng)
    point = setprecision(() -> _testhelper_consistent_point(m, rng), BigFloat, 256)
    concs = Dict{Symbol, BigFloat}(x => exp(3 * big(rand(rng)) - big(1.5))
                                   for x in EnzymeRates.metabolites(m))
    params = Dict{Symbol, BigFloat}(k => point.vals[k] for k in fitted)
    params[:Keq] = point.Keq
    point, concs, params
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
    # The point and the solver take the concrete mechanism: each compiled mechanism is
    # its own type, so given one they would compile again for every mechanism.
    mech = _testhelper_lift(em)
    fitted = EnzymeRates._dependent_param_exprs(typeof(em))[2]
    mets = EnzymeRates.metabolites(em)
    for _ in 1:n_points
        point, concs, params = _testhelper_law_inputs(mech, fitted, rng)
        oracle = _testhelper_steady_state_rate(mech, point, concs)
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

# ── Well-posedness oracles ──────────────────────────────────────────────────
#
# The well-posedness guard's oracles (spec section 14): mass-action networks of one
# catalytic subunit, read from the graph read API and never from the guard's code, whose
# conformations flip at free E (form 1); the RE-speed oracle on the guard's faces; the
# ε-path oracle of admission, along integer potentials or the corners of the admissible
# cone and five orders of the two limits; the exact LP of admission's existence test; and
# the concerted n-protomer oracle.

"""
    @_testhelper_unguarded call

The concrete mechanism that `call`, an `@enzyme_mechanism` or `@allosteric_mechanism`
call, writes: the constructors on the DSL walker's raw value, without the well-posedness
check the macros run, so a test can hold a mechanism the macros refuse.
"""
macro _testhelper_unguarded(call)
    allosteric, raw = _testhelper_mechanism_call(call, "@_testhelper_unguarded")
    esc(allosteric ? :(EnzymeRates.AllostericMechanism($raw)) :
                     :(EnzymeRates.Mechanism($raw)))
end

"""
The faces of `m` (concrete or compiled) that the guard checks (spec section 8), each as
the set of names it zeroes, in the guard's order: the interior, every product zero, every
substrate zero, each substrate and then each product alone on a side of two or more; then
each of these with every modifier-only name zero too. A reactant's `::Mod` copy carries
its reactant's name, so it is zeroed with it; a regulatory-site ligand never is.
"""
function _testhelper_faces(m)
    cm = EnzymeRates.catalytic(_testhelper_lift(m))
    mods = Set(x.name for x in cm.table
               if x.modifier && !(x.name in cm.substrates || x.name in cm.products))
    base = [Set{Symbol}(), Set(cm.products), Set(cm.substrates)]
    length(cm.substrates) >= 2 && append!(base, [Set([x]) for x in cm.substrates])
    length(cm.products) >= 2 && append!(base, [Set([x]) for x in cm.products])
    isempty(mods) ? base : [base; [union(z, mods) for z in base]]
end

"""Whether conformation `c` (1 active, 2 inactive) of `m` excludes edge `k` of its
catalytic graph: an `:OnlyA` group's edge in the inactive conformation, an `:OnlyI`
group's in the active one. A plain mechanism excludes nothing."""
function _testhelper_excludes(m, k::Int, c::Int)
    m isa EnzymeRates.AllostericMechanism || return false
    tag = EnzymeRates.cat_allo_state(m, m.catalytic.group[k])
    (c == 2 && tag === :OnlyA) || (c == 1 && tag === :OnlyI)
end

"""One reversible step of an oracle network between nodes `i` and `j`: forward rate `f`
and reverse rate `r`, concentrations included, `w` copies of the first substrate taken up
from `i` to `j`; a `fast` link runs at Λ times its rates."""
struct _testhelper_Link
    i::Int
    j::Int
    f::BigFloat
    r::BigFloat
    w::Int
    fast::Bool
end

"""A mass-action network of the oracles: its nodes, each `(form, conformation)` (`(state,
conformation)` for the n-protomer oracle's product states), and its links."""
struct _testhelper_Network
    nodes::Vector{Tuple{Int, Int}}
    links::Vector{_testhelper_Link}
end

"""The equilibrium ratio of the flip at free E, inactive over active, of the allosteric
mechanism `mech` at `concs`: `L` times the ratio of the two conformations'
regulatory-site factors. Regulators bind each conformation at rapid equilibrium whatever
its catalytic form, so summed over regulator occupancy the flip takes this ratio."""
function _testhelper_flip_ratio(mech, point, concs::AbstractDict)
    function _testhelper_site_factor(inactive)
        prod((1 + sum((concs[lig] / point.vals[K]
                       for (lig, tag) in zip(site.ligands, site.allo_states)
                       for K in (EnzymeRates._reg_K(mech, lig, tag, inactive),)
                       if K !== nothing); init = big(0)))^site.multiplicity
             for site in EnzymeRates.regulatory_sites(mech); init = big(1))
    end
    point.L * _testhelper_site_factor(true) / _testhelper_site_factor(false)
end

"""Speed factors of the rapid-equilibrium edges, one per edge of `m`'s catalytic graph and
conformation, log-uniform over `10^±3`: a speed pattern of the RE-speed oracle."""
_testhelper_draw_speeds(m, rng) =
    BigFloat[big(10)^(3 * (2 * big(rand(rng)) - 1))
             for _ in EnzymeRates.catalytic(_testhelper_lift(m)).edges, _ in 1:2]

"""
    _testhelper_network(m, point, concs; speeds, mode = :law, ε = 1, potentials, sigma,
                        theta) -> _testhelper_Network

The mass-action network of `m`'s catalytic graph, one catalytic subunit, at the constants
of `point` (`_testhelper_consistent_point`) and the concentrations `concs`: one copy of
the forms for a plain mechanism, two for an allosteric one, joined by a fast flip between
the copies of free E at `_testhelper_flip_ratio`. A rapid-equilibrium edge of
conformation `c` runs fast at speed `speeds[k, c]`, a steady-state one at
`point.kf[k, c]`.

`mode = :law` is the law's network: each conformation keeps the edges it does not
exclude, on the forms they join to free E (deletion and pruning). `mode = :vanish` is the
ε-network: every form and edge, each excluded edge vanishing as ε → 0 along a direction,
with `potentials[c][f]` the potential of form `f`'s component in conformation `c` and
`d = potentials[c][to] - potentials[c][from]`. A plain binding's equilibrium constant
scales by `ε^d`: a rapid-equilibrium one stays fast, its on rate scaled by `ε^(θ d)` and
its off rate by `ε^(-(1 - θ) d)`; a steady-state one keeps its off rate. Any other
excluded edge is slow, its rates scaled by `ε^(σ + max(±d, 0))`. `sigma` (σ) and `theta`
(θ) are one value for every edge or one per edge.
"""
function _testhelper_network(m, point, concs::AbstractDict; speeds::Matrix{BigFloat},
                             mode::Symbol = :law, ε::BigFloat = big(1.0),
                             potentials::Vector{Vector{BigFloat}} = Vector{BigFloat}[],
                             sigma::Union{BigFloat, Vector{BigFloat}} = big(1.0),
                             theta::Union{BigFloat, Vector{BigFloat}} = big(1.0) / 2)
    mech = _testhelper_lift(m)
    cm = EnzymeRates.catalytic(mech)
    nf = length(cm.forms)
    nc = mech isa EnzymeRates.AllostericMechanism ? 2 : 1
    x1 = first(cm.substrates)
    met_names(c) = [cm.table[i].name for i in EnzymeRates.entries(c)]
    keep = [trues(nf) for _ in 1:nc]
    if mode === :law
        for c in 1:nc
            reached = falses(nf)
            reached[1] = true
            grew = true
            while grew
                grew = false
                for (k, e) in enumerate(cm.edges)
                    !_testhelper_excludes(mech, k, c) &&
                        reached[e.from] != reached[e.to] || continue
                    reached[e.from] = reached[e.to] = true
                    grew = true
                end
            end
            keep[c] = reached
        end
    end
    nodes = [(f, c) for c in 1:nc for f in 1:nf if keep[c][f]]
    id = Dict(node => n for (n, node) in enumerate(nodes))
    links = _testhelper_Link[]
    for c in 1:nc, (k, e) in enumerate(cm.edges)
        keep[c][e.from] && keep[c][e.to] || continue
        excluded = _testhelper_excludes(mech, k, c)
        mode === :law && excluded && continue
        K, kf = point.K[k, c], point.kf[k, c]
        fast = e.rapid_equilibrium
        f, r = fast ? (speeds[k, c] * sqrt(K), speeds[k, c] / sqrt(K)) : (kf, kf / K)
        if excluded
            d = potentials[c][e.to] - potentials[c][e.from]
            if EnzymeRates._is_plain_binding(cm, e)
                θ = theta isa BigFloat ? theta : theta[k]
                f, r = fast ? (f * ε^(θ * d), r * ε^(-(1 - θ) * d)) : (f * ε^d, r)
            else
                σ = sigma isa BigFloat ? sigma : sigma[k]
                f, r, fast = f * ε^(σ + max(d, 0)), r * ε^(σ + max(-d, 0)), false
            end
        end
        w = count(==(x1), met_names(e.consumed)) - count(==(x1), met_names(e.released))
        conc_in = prod((concs[x] for x in met_names(e.consumed)); init = big(1))
        conc_out = prod((concs[x] for x in met_names(e.released)); init = big(1))
        push!(links, _testhelper_Link(id[(e.from, c)], id[(e.to, c)], f * conc_in,
                                      r * conc_out, w, fast))
    end
    if nc == 2
        ratio = _testhelper_flip_ratio(mech, point, concs)
        push!(links, _testhelper_Link(id[(1, 1)], id[(1, 2)], sqrt(ratio),
                                      1 / sqrt(ratio), 0, true))
    end
    _testhelper_Network(nodes, links)
end

"""The closed classes (strongly connected, no transition leaving) of the digraph on nodes
`1:n` with transitions `arcs`, by Tarjan's algorithm, each as its sorted nodes."""
function _testhelper_closed_classes(n::Int, arcs::Vector{Tuple{Int, Int}})
    next = [Int[] for _ in 1:n]
    foreach(((a, b),) -> push!(next[a], b), arcs)
    index, low, comp = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    onstack, stack = falses(n), Int[]
    counter, ncomp = 0, 0
    for s in 1:n
        index[s] == 0 || continue
        counter += 1
        index[s] = low[s] = counter
        push!(stack, s)
        onstack[s] = true
        work = [(s, 1)]
        while !isempty(work)
            u, i = work[end]
            if i <= length(next[u])
                work[end] = (u, i + 1)
                v = next[u][i]
                if index[v] == 0
                    counter += 1
                    index[v] = low[v] = counter
                    push!(stack, v)
                    onstack[v] = true
                    push!(work, (v, 1))
                elseif onstack[v]
                    low[u] = min(low[u], index[v])
                end
            else
                pop!(work)
                isempty(work) || (low[work[end][1]] = min(low[work[end][1]], low[u]))
                if low[u] == index[u]
                    ncomp += 1
                    while true
                        v = pop!(stack)
                        onstack[v] = false
                        comp[v] = ncomp
                        v == u && break
                    end
                end
            end
        end
    end
    leaves = falses(ncomp)
    for (a, b) in arcs
        comp[a] != comp[b] && (leaves[comp[a]] = true)
    end
    [findall(==(k), comp) for k in 1:ncomp if !leaves[k]]
end

"""The stationary distribution of the irreducible chain with off-diagonal rates `R[i, j]`
(`i → j`), by Grassmann–Taksar–Heyman state reduction: sums of non-negative terms only, so
every entry keeps its relative precision whatever the spread of the rates."""
function _testhelper_gth(R::Matrix{BigFloat})
    n = size(R, 1)
    A = copy(R)
    for k in n:-1:2
        s = sum(A[k, j] for j in 1:(k - 1))
        iszero(s) && error("_testhelper_gth: state $k has no exit to states 1:$(k - 1)")
        for i in 1:(k - 1)
            iszero(A[i, k]) && continue
            A[i, k] /= s
            for j in 1:(k - 1)
                i == j || iszero(A[k, j]) || (A[i, j] += A[i, k] * A[k, j])
            end
        end
    end
    p = zeros(BigFloat, n)
    p[1] = 1
    for j in 2:n
        p[j] = sum((p[i] * A[i, j] for i in 1:(j - 1) if !iszero(A[i, j])); init = big(0))
    end
    p ./ sum(p)
end

"""Each node's potential of the first substrate along `links` (node `j`'s is node `i`'s
plus the link's uptake `w`), component by component on nodes `1:n`; `nothing` when a
cycle of the links takes the substrate up."""
function _testhelper_offsets(n::Int, links)
    h = fill(typemin(Int), n)
    adjacent = [Tuple{Int, Int}[] for _ in 1:n]
    for l in links
        push!(adjacent[l.i], (l.j, l.w))
        push!(adjacent[l.j], (l.i, -l.w))
    end
    for s in 1:n
        h[s] == typemin(Int) || continue
        h[s] = 0
        todo = [s]
        while !isempty(todo)
            u = pop!(todo)
            for (v, w) in adjacent[u]
                if h[v] == typemin(Int)
                    h[v] = h[u] + w
                    push!(todo, v)
                elseif h[v] != h[u] + w
                    return nothing
                end
            end
        end
    end
    h
end

"""
    _testhelper_steady_rate(nw; Λ) -> (; v, nclasses)

The steady-state turnover per enzyme of `nw` (net uptake of the first substrate) with
every fast link at `Λ` times its rates, in the current BigFloat precision; `v` is NaN when
the chain has more than one closed class (`nclasses` of them). It weights each link's flux
by its uptake plus the fall of the first substrate's potential along the fast links that
run both ways, so those links, whose fluxes nearly cancel at large Λ, drop out; when they
close a cycle that takes the substrate up, every link counts with its own uptake.
"""
function _testhelper_steady_rate(nw::_testhelper_Network; Λ::BigFloat)
    n = length(nw.nodes)
    rates(l) = l.fast ? (Λ * l.f, Λ * l.r) : (l.f, l.r)
    arcs = Tuple{Int, Int}[]
    for l in nw.links
        l.f > 0 && push!(arcs, (l.i, l.j))
        l.r > 0 && push!(arcs, (l.j, l.i))
    end
    classes = _testhelper_closed_classes(n, arcs)
    length(classes) == 1 || return (v = big(NaN), nclasses = length(classes))
    C = only(classes)
    at = zeros(Int, n)
    at[C] .= eachindex(C)
    inside = [l for l in nw.links if at[l.i] > 0 && at[l.j] > 0]
    R = zeros(BigFloat, length(C), length(C))
    for l in inside
        rf, rr = rates(l)
        R[at[l.i], at[l.j]] += rf
        R[at[l.j], at[l.i]] += rr
    end
    p = length(C) == 1 ? BigFloat[1] : _testhelper_gth(R)
    two_way(l) = l.fast && l.f > 0 && l.r > 0
    h = _testhelper_offsets(n, filter(two_way, inside))
    v = big(0.0)
    for l in inside
        rf, rr = rates(l)
        k = h === nothing ? l.w : two_way(l) ? 0 : l.w + h[l.i] - h[l.j]
        k == 0 || (v += k * (rf * p[at[l.i]] - rr * p[at[l.j]]))
    end
    (v = v, nclasses = 1)
end

"""
The rapid-equilibrium limit (Λ = ∞) of `nw`, exactly, at an interior point: the forms of
each segment of fast links take their equilibrium weights, the segments form the chain of
the slow links, and the turnover sums the slow links as `_testhelper_steady_rate` does.
`Inf` when the fast links close a cycle that takes the first substrate up, NaN when the
slow chain has more than one closed class.
"""
function _testhelper_re_limit(nw::_testhelper_Network)
    n = length(nw.nodes)
    fast = [l for l in nw.links if l.fast]
    all(l -> l.f > 0 && l.r > 0, fast) ||
        error("_testhelper_re_limit: a fast link has a zero rate; the limit is taken at " *
              "interior points")
    h = _testhelper_offsets(n, fast)
    h === nothing && return big(Inf)
    segment, weight = zeros(Int, n), zeros(BigFloat, n)
    adjacent = [Tuple{Int, BigFloat}[] for _ in 1:n]
    for l in fast
        push!(adjacent[l.i], (l.j, l.f / l.r))
        push!(adjacent[l.j], (l.i, l.r / l.f))
    end
    ns = 0
    for s in 1:n
        segment[s] == 0 || continue
        ns += 1
        segment[s], weight[s] = ns, big(1)
        todo = [s]
        while !isempty(todo)
            u = pop!(todo)
            for (v, q) in adjacent[u]
                segment[v] == 0 || continue
                segment[v], weight[v] = ns, weight[u] * q
                push!(todo, v)
            end
        end
    end
    W = zeros(BigFloat, ns)
    foreach(i -> W[segment[i]] += weight[i], 1:n)
    R, arcs = zeros(BigFloat, ns, ns), Tuple{Int, Int}[]
    for l in nw.links
        l.fast && continue
        a, b = segment[l.i], segment[l.j]
        a == b && continue
        rf, rr = l.f * weight[l.i] / W[a], l.r * weight[l.j] / W[b]
        rf > 0 && (R[a, b] += rf; push!(arcs, (a, b)))
        rr > 0 && (R[b, a] += rr; push!(arcs, (b, a)))
    end
    classes = _testhelper_closed_classes(ns, arcs)
    length(classes) == 1 || return big(NaN)
    C = only(classes)
    P = zeros(BigFloat, ns)
    P[C] .= length(C) == 1 ? BigFloat[1] : _testhelper_gth(R[C, C])
    p = [P[segment[i]] * weight[i] / W[segment[i]] for i in 1:n]
    v = big(0.0)
    for l in nw.links
        l.fast && continue
        k = l.w + h[l.i] - h[l.j]
        k == 0 || (v += k * (l.f * p[l.i] - l.r * p[l.j]))
    end
    v
end

"""The relative gap between `x` and `y`, with an absolute floor `least`."""
_testhelper_relgap(x, y; least = big(10)^-14) = abs(x - y) / max(abs(x), abs(y), least)

"""The package's rate law of the concrete mechanism `mech` at the constants of `point`
(its fitted values and `Keq`) and the concentrations `concs`, read without the boundary
check (`_testhelper_law_value`); NaN where it has no value."""
function _testhelper_package_law(mech, point, concs::AbstractDict)
    params = Dict{Symbol, BigFloat}(k => point.vals[k]
                                    for k in EnzymeRates._dependent_param_exprs(mech)[2])
    params[:Keq] = point.Keq
    _testhelper_law_value(mech, params, concs)
end

"""
The law's value on the face of `mech` where the names `face` vanish, computed without the
package: the rapid-equilibrium limit of the law's own network (`_testhelper_re_limit`) at
the interior points `face = δ·u`, δ = 1e-40, 1e-80, 1e-160, along two random directions
`u`; NaN when it grows without bound, has not settled by δ = 1e-160 or depends on `u` (the
law's 0/0 does not cancel). With `face` empty, the limit at `concs`.
"""
function _testhelper_limit_law(mech, point, concs::AbstractDict, face::Set{Symbol}; rng)
    speeds = _testhelper_draw_speeds(mech, rng)
    value(c) = _testhelper_re_limit(_testhelper_network(mech, point, c; speeds))
    inner = value(concs)
    isempty(face) && return inner
    tiny = big(10)^-30 * max(abs(inner), big(10)^-40)
    limits = BigFloat[]
    for _ in 1:2
        u = Dict(z => big(10)^(2 * big(rand(rng)) - 1) for z in face)
        vals = [value(merge(concs, Dict(z => δ * u[z] for z in face)))
                for δ in (big(10)^-40, big(10)^-80, big(10)^-160)]
        all(isfinite, vals) || return big(NaN)
        abs(vals[3]) > 1e10 * max(abs(vals[2]), tiny) && return big(NaN)
        _testhelper_relgap(vals[2], vals[3]; least = tiny) > big(10)^-25 && return big(NaN)
        push!(limits, vals[3])
    end
    _testhelper_relgap(limits[1], limits[2]; least = tiny) > big(10)^-20 ? big(NaN) :
        limits[1]
end

"""
    _testhelper_respeed(m, point, concs, face; rng, law = :package) -> Symbol

The RE-speed oracle (spec section 14) on the face of `m` where the names `face` vanish:
the law's network (`_testhelper_network`, one subunit) solved at Λ = 1e10 and 1e20 under
two speed patterns of 10^±3, its rapid-equilibrium edges at Λ times their pattern's
factor. An allosteric `m` is judged at one catalytic subunit (`_testhelper_single_subunit`).
The first verdict that applies:
- `:nonunique`: the face chain has more than one closed class;
- `:diverges`: |v| grows more than 1e3-fold from the first Λ to the last (RE turnover);
- `:unconverged`: the change from the first Λ to the last, scaled by their ratio, exceeds
  1e-10 relative, so the last Λ is too small to judge 1e-8;
- `:pattern`: the values at the last Λ differ between patterns by more than 1e-8;
- `:law_nan`: the law has no value on the face;
- `:mismatch`: the limit differs from the law by more than 1e-8 relative;
- `:agree`.
The law is the package's (`law = :package`, read unchecked) or the oracle's own
(`law = :limit`, `_testhelper_limit_law`).
"""
function _testhelper_respeed(m, point, concs::AbstractDict, face::Set{Symbol}; rng,
                             law::Symbol = :package)
    mech = _testhelper_single_subunit(_testhelper_lift(m))
    Λs = (big(10.0)^10, big(10.0)^20)
    c0 = merge(concs, Dict(z => big(0) for z in face))
    patterns = [_testhelper_draw_speeds(mech, rng) for _ in 1:2]
    solved = [_testhelper_steady_rate(_testhelper_network(mech, point, c0; speeds); Λ)
              for speeds in patterns, Λ in Λs]
    any(x -> x.nclasses > 1, solved) && return :nonunique
    v = [x.v for x in solved]
    top = size(v, 2)
    any(i -> abs(v[i, top]) > 1e3 * max(abs(v[i, 1]), big(10)^-30), axes(v, 1)) &&
        return :diverges
    any(i -> abs(v[i, 1] - v[i, top]) * Λs[1] / Λs[end] >
             1e-10 * max(abs(v[i, top]), big(10)^-14), axes(v, 1)) && return :unconverged
    any(i -> _testhelper_relgap(v[i, top], v[1, top]) > 1e-8, axes(v, 1)) &&
        return :pattern
    vlaw = law === :package ? _testhelper_package_law(mech, point, c0) :
           _testhelper_limit_law(mech, point, concs, face; rng)
    isfinite(vlaw) || return :law_nan
    _testhelper_relgap(v[1, top], vlaw) > 1e-8 ? :mismatch : :agree
end

"""
The RE-speed oracle's verdict on each face of `m` (`_testhelper_faces`), the worst over
`npoints` points (`_testhelper_consistent_point`, concentrations log-uniform over
10^±1) in the order nonunique, diverges, unconverged, pattern, law_nan, mismatch, agree.
An allosteric mechanism is judged at one catalytic subunit (`_testhelper_single_subunit`).
"""
function _testhelper_respeed_verdicts(m; rng, npoints::Int = 1, law::Symbol = :package)
    mech = _testhelper_single_subunit(_testhelper_lift(m))
    faces = _testhelper_faces(mech)
    order = (:nonunique, :diverges, :unconverged, :pattern, :law_nan, :mismatch, :agree)
    worst = fill(:agree, length(faces))
    setprecision(BigFloat, 256) do
        for _ in 1:npoints
            point = _testhelper_consistent_point(mech, rng)
            concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                           for x in EnzymeRates.metabolites(mech))
            for (i, face) in enumerate(faces)
                v = _testhelper_respeed(mech, point, concs, face; rng, law)
                findfirst(==(v), order) < findfirst(==(worst[i]), order) &&
                    (worst[i] = v)
            end
        end
    end
    worst
end

"""
    _testhelper_exclusion(am, c) -> (; component, arcs, failure, lost)

The exclusion structure of conformation `c` of `am` (spec section 7): `component[f]`, the
component of form `f` under the edges `c` keeps, numbered by first form; `arcs`,
`(unbound side's component, bound side's, edge)` for each plain binding `c` excludes;
`failure`, `:inside` when an arc stays inside one component, `:cycle` when the arcs close
a directed cycle, else `:none` (a path to the limit exists); `lost`, the components free
E's does not reach along the arcs (islands).
"""
function _testhelper_exclusion(am, c::Int)
    cm = EnzymeRates.catalytic(am)
    nf = length(cm.forms)
    parent = collect(1:nf)
    root(i) = parent[i] == i ? i : root(parent[i])
    for (k, e) in enumerate(cm.edges)
        _testhelper_excludes(am, k, c) || (parent[root(e.from)] = root(e.to))
    end
    roots = [root(f) for f in 1:nf]
    label = Dict(r => q for (q, r) in enumerate(unique(roots)))
    component = [label[r] for r in roots]
    n = length(label)
    arcs = Tuple{Int, Int, Int}[(component[e.from], component[e.to], k)
                                for (k, e) in enumerate(cm.edges)
                                if _testhelper_excludes(am, k, c) &&
                                   EnzymeRates._is_plain_binding(cm, e)]
    failure = any(a -> a[1] == a[2], arcs) ? :inside : :none
    if failure === :none
        indegree = zeros(Int, n)
        foreach(a -> indegree[a[2]] += 1, arcs)
        todo = findall(iszero, indegree)
        while !isempty(todo)
            u = pop!(todo)
            for a in arcs
                a[1] == u && (indegree[a[2]] -= 1) == 0 && push!(todo, a[2])
            end
        end
        any(>(0), indegree) && (failure = :cycle)
    end
    reached = falses(n)
    reached[component[1]] = true
    grew = true
    while grew
        grew = false
        for a in arcs
            reached[a[1]] && !reached[a[2]] && (reached[a[2]] = true; grew = true)
        end
    end
    (; component, arcs, failure, lost = findall(!, reached))
end

"""The excluding conformations of `am` with their exclusion structures, `[(c, ex)]`."""
function _testhelper_excluding(am)
    edges = eachindex(EnzymeRates.catalytic(am).edges)
    [(c, _testhelper_exclusion(am, c)) for c in 1:2
     if any(k -> _testhelper_excludes(am, k, c), edges)]
end

"""Potentials of the components of `ex` (`_testhelper_exclusion` with no failure): each
starts at `start[x]`, each arc raises its bound side to at least its unbound side's plus
`rise[q]`, and free E's component is shifted back to 0."""
function _testhelper_component_potentials(ex, start::Vector{BigFloat},
                                          rise::Vector{BigFloat})
    φ = copy(start)
    for _ in eachindex(φ), (q, a) in enumerate(ex.arcs)
        φ[a[2]] = max(φ[a[2]], φ[a[1]] + rise[q])
    end
    φ .- φ[ex.component[1]]
end

"""
The ε-path's six standard directions for `am`: integer potentials, arcs rising by 1 and
every component other than free E's starting at +1 or at -1, each with one speed exponent
σ ∈ {½, 1, 2} for every blocked edge and the even split θ = ½ of each vanishing
rapid-equilibrium binding. Each is `(; potentials, sigma, theta)`, `potentials[c][f]` the
potential of form `f` in conformation `c` (empty where `c` excludes nothing).
"""
function _testhelper_standard_directions(am)
    ne = length(EnzymeRates.catalytic(am).edges)
    directions = []
    for offset in (1, -1), σ in (big(0.5), big(1.0), big(2.0))
        potentials = [BigFloat[] for _ in 1:2]
        for (c, ex) in _testhelper_excluding(am)
            start = [x == ex.component[1] ? big(0.0) : BigFloat(offset)
                     for x in 1:maximum(ex.component)]
            φ = _testhelper_component_potentials(ex, start, ones(BigFloat, length(ex.arcs)))
            potentials[c] = [φ[x] for x in ex.component]
        end
        push!(directions, (; potentials, sigma = fill(σ, ne), theta = fill(big(0.5), ne)))
    end
    directions
end

"""The corners of a box of `k` binary choices, each a vector of `k` zeros and ones: all
`2^k` when there are at most 128, else 128 distinct ones drawn from `rng`."""
function _testhelper_corners(k::Int, rng)
    k == 0 && return [Int[]]
    k <= 7 && return [digits(x; base = 2, pad = k) for x in 0:(2^k - 1)]
    drawn = Set{Vector{Int}}()
    while length(drawn) < 128
        push!(drawn, rand(rng, 0:1, k))
    end
    sort!(collect(drawn))
end

"""
The ε-path's corner directions for `am`: the corners of the admissible cone, every binary
choice made one way or the other. Every arc rises by 1 or 20, every component other than
free E's that no arc enters starts at -5 or 5, every blocked edge's σ is 1/8 or 16 and
every vanishing rapid-equilibrium binding's θ is 0 or 1; all corners when there are at
most 128, else 128 of them (`_testhelper_corners`).
"""
function _testhelper_corner_directions(am, rng)
    cm = EnzymeRates.catalytic(am)
    excluding = _testhelper_excluding(am)
    slots = Tuple{Symbol, Int, Int}[]
    for (c, ex) in excluding
        append!(slots, [(:rise, c, q) for q in eachindex(ex.arcs)])
        entered = Set(a[2] for a in ex.arcs)
        append!(slots, [(:start, c, x) for x in 1:maximum(ex.component)
                        if x != ex.component[1] && !(x in entered)])
    end
    for (k, e) in enumerate(cm.edges)
        any(c -> _testhelper_excludes(am, k, c), 1:2) || continue
        if !EnzymeRates._is_plain_binding(cm, e)
            push!(slots, (:sigma, k, 0))
        elseif e.rapid_equilibrium
            push!(slots, (:theta, k, 0))
        end
    end
    map(_testhelper_corners(length(slots), rng)) do bits
        high(slot) = (q = findfirst(==(slot), slots); q !== nothing && bits[q] == 1)
        potentials = [BigFloat[] for _ in 1:2]
        for (c, ex) in excluding
            entered = Set(a[2] for a in ex.arcs)
            start = [x == ex.component[1] ? big(0.0) : x in entered ? big(-Inf) :
                     high((:start, c, x)) ? big(5.0) : big(-5.0)
                     for x in 1:maximum(ex.component)]
            rise = [high((:rise, c, q)) ? big(20.0) : big(1.0) for q in eachindex(ex.arcs)]
            φ = _testhelper_component_potentials(ex, start, rise)
            potentials[c] = [φ[x] for x in ex.component]
        end
        edges = eachindex(cm.edges)
        (; potentials,
         sigma = [high((:sigma, k, 0)) ? big(16.0) : big(0.125) for k in edges],
         theta = [high((:theta, k, 0)) ? big(1.0) : big(0.0) for k in edges])
    end
end

"""How far the ε-network of `am` along direction `d` spreads its rates' ε-exponents
(`_testhelper_network`): the sum over the excluded edges of the magnitudes of each edge's
two rate exponents. With δ the potential's rise across the edge and σ its speed exponent,
that is |δ| for a plain binding and 2σ + |δ| for any other edge, whose exponents are
σ + max(δ, 0) and σ + max(-δ, 0). The RE-first order lets Λ outgrow it."""
function _testhelper_epsilon_span(am, d)
    cm = EnzymeRates.catalytic(am)
    span = big(0.0)
    for c in 1:2, (k, e) in enumerate(cm.edges)
        _testhelper_excludes(am, k, c) || continue
        rise = abs(d.potentials[c][e.to] - d.potentials[c][e.from])
        span += EnzymeRates._is_plain_binding(cm, e) ? rise : 2 * d.sigma[k] + rise
    end
    span
end

"""
The limit of the turnover of the ε-network `network(ε)` along Λ = ε^-κ as ε → 0, ε =
10^-E with E doubling: `(; v, converged)`. κ is `:re_first` (Λ outgrows every ε-scale,
`2 span + 10`), `:eps_first` (Λ grows slower than any rate vanishes, `1/(8(nodes + 1))`)
or a number. Converged once Λ ≥ 1e40 and the last two steps move by less than 1e-12 and
1e-6 of `scale`; infinite when the rate keeps growing by orders of magnitude.
"""
function _testhelper_path_limit(network, κ, span::BigFloat, nnodes::Int, scale::BigFloat)
    k = κ === :re_first ? 2 * span + 10 :
        κ === :eps_first ? big(1) / (8 * (nnodes + 1)) : BigFloat(κ)
    E, Emax = κ === :eps_first ? (4, 2^18) : (2, 2^14)
    vals = BigFloat[]
    while E <= Emax
        ε = big(10.0)^(-E)
        push!(vals, _testhelper_steady_rate(network(ε); Λ = ε^(-k)).v)
        n = length(vals)
        if n >= 3
            d1, d2 = abs(vals[n] - vals[n - 1]), abs(vals[n - 1] - vals[n - 2])
            k * E >= 40 && d1 <= big(10)^-12 * scale && d2 <= big(10)^-6 * scale &&
                return (v = vals[n], converged = true)
            abs(vals[n]) > 1e6 * max(abs(vals[n - 1]), scale) &&
                abs(vals[n - 1]) > 1e3 * max(abs(vals[n - 2]), scale) &&
                return (v = big(Inf), converged = true)
        end
        E *= 2
    end
    (v = vals[end], converged = false)
end

"""
    _testhelper_epath_verdict(m; rng, face = Set{Symbol}(), corners = false,
                              precision = 256) -> Symbol

The ε-path oracle (spec section 14) on the face of `m` where the names `face` vanish, one
catalytic subunit: along each direction (the six standard ones, or the corners with
`corners`) and each order κ ∈ (RE first, 8, 2, ½, ε first) of the two limits, the
turnover of the ε-network against the deletion law's value there
(`_testhelper_limit_law`). The first verdict that applies:
- `:no_path`: an excluding conformation has no path to the limit (`_testhelper_exclusion`);
- `:nonunique`: an ε-network has more than one closed class on the face;
- `:diverges`: some path's turnover grows without bound (a fast cycle of the ε-network
  turns over);
- `:del_nonunique`: the deletion law has no value on the face;
- `:converges`: every path tends to it, to 1e-8 of the rate's scale;
- `:wrong_limit`: every path settles, none on it, all on one other value;
- `:direction_dependent`: otherwise (some paths tend to it and some do not, or some never
  settle).
"""
function _testhelper_epath_verdict(m; rng, face::Set{Symbol} = Set{Symbol}(),
                                   corners::Bool = false, precision::Int = 256)
    setprecision(BigFloat, precision) do
        am = _testhelper_single_subunit(_testhelper_lift(m))
        any(((_, ex),) -> ex.failure !== :none, _testhelper_excluding(am)) &&
            return :no_path
        point = _testhelper_consistent_point(am, rng)
        concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(am))
        c0 = merge(concs, Dict(z => big(0) for z in face))
        vdel = _testhelper_limit_law(am, point, concs, face; rng)
        vin = _testhelper_limit_law(am, point, concs, Set{Symbol}(); rng)
        scale = max(isfinite(vdel) ? abs(vdel) : big(0), isfinite(vin) ? abs(vin) : big(1),
                    big(10)^-30)
        directions = corners ? _testhelper_corner_directions(am, rng) :
                     _testhelper_standard_directions(am)
        speeds = _testhelper_draw_speeds(am, rng)
        outcomes, ends = Symbol[], BigFloat[]
        for d in directions
            function _testhelper_eps_network(ε)
                _testhelper_network(am, point, c0; speeds, mode = :vanish, ε,
                                    potentials = d.potentials, sigma = d.sigma,
                                    theta = d.theta)
            end
            nw = _testhelper_eps_network(big(10)^-2)
            _testhelper_steady_rate(nw; Λ = big(1.0)).nclasses > 1 && return :nonunique
            span = _testhelper_epsilon_span(am, d)
            for κ in (:re_first, 8, 2, 1 // 2, :eps_first)
                r = _testhelper_path_limit(_testhelper_eps_network, κ, span,
                                           length(nw.nodes), scale)
                if !r.converged
                    push!(outcomes, :unsettled)
                elseif !isfinite(r.v)
                    push!(outcomes, :infinite)
                elseif abs(r.v - vdel) <= big(10)^-8 * scale
                    push!(outcomes, :tends)
                else
                    push!(outcomes, :elsewhere)
                    push!(ends, r.v)
                end
            end
        end
        :infinite in outcomes && return :diverges
        isfinite(vdel) || return :del_nonunique
        all(==(:tends), outcomes) && return :converges
        all(==(:elsewhere), outcomes) &&
            all(x -> abs(x - ends[1]) <= big(10)^-6 * scale, ends) && return :wrong_limit
        :direction_dependent
    end
end

"""
Whether every standard direction (`_testhelper_standard_directions`) of `m`'s ε-path at an
interior point approaches the deletion law at a power of ε: over ε = 1e-9 … 1e-18 the
error relative to the rate is `C·ε^p` within 10% at every decade, `p` a positive multiple
of ½, with exact rapid-equilibrium limits (`_testhelper_re_limit`) in 512-bit BigFloat; an
error below 1e-60 throughout passes.
"""
function _testhelper_in_epath_window(m; rng)
    setprecision(BigFloat, 512) do
        am = _testhelper_single_subunit(_testhelper_lift(m))
        point = _testhelper_consistent_point(am, rng)
        concs = Dict{Symbol, BigFloat}(x => big(10)^(2 * big(rand(rng)) - 1)
                                       for x in EnzymeRates.metabolites(am))
        speeds = _testhelper_draw_speeds(am, rng)
        vdel = _testhelper_re_limit(_testhelper_network(am, point, concs; speeds))
        scale = max(abs(vdel), big(10)^-30)
        eps = [big(10.0)^-j for j in 9:18]
        all(_testhelper_standard_directions(am)) do d
            errors = [abs(_testhelper_re_limit(_testhelper_network(
                          am, point, concs; speeds, mode = :vanish, ε,
                          potentials = d.potentials, sigma = d.sigma,
                          theta = d.theta)) - vdel) / scale for ε in eps]
            all(<(big(10)^-60), errors) && return true
            any(iszero, errors) && return false
            slopes = sort([log10(errors[j] / errors[j + 1]) for j in 1:(length(eps) - 1)])
            p = slopes[(length(slopes) + 1) ÷ 2]
            q = round(2 * p) / 2
            q >= 1 // 2 && abs(p - q) < 0.05 || return false
            C = errors[end] / eps[end]^q
            all(j -> abs(errors[j] / (C * eps[j]^q) - 1) < 0.1, eachindex(eps))
        end
    end
end

"""Whether `{x ≥ 0 : M x = b}` is nonempty: phase 1 of the simplex method with Bland's
rule, in exact rational arithmetic."""
function _testhelper_feasible(M::Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    m, n = size(M)
    flip = [x < 0 ? -1 : 1 for x in b]
    T = hcat(flip .* M, Matrix{Rational{BigInt}}(LinearAlgebra.I, m, m), flip .* b)
    basis = collect((n + 1):(n + m))
    cost = [zeros(Rational{BigInt}, n); ones(Rational{BigInt}, m)]
    while true
        reduced(j) = cost[j] - sum(cost[basis[i]] * T[i, j] for i in 1:m)
        j = findfirst(j -> reduced(j) < 0, 1:(n + m))
        j === nothing && break
        rows = [i for i in 1:m if T[i, j] > 0]
        ratio = minimum(i -> T[i, end] / T[i, j], rows)
        i = argmin(i -> basis[i], [i for i in rows if T[i, end] / T[i, j] == ratio])
        T[i, :] ./= T[i, j]
        for k in 1:m
            k == i || iszero(T[k, j]) || (T[k, :] .-= T[k, j] .* T[i, :])
        end
        basis[i] = j
    end
    iszero(sum(cost[basis[i]] * T[i, end] for i in 1:m))
end

"""
Admission's existence test as an exact LP (spec section 7): whether conformation `c` of
`am` has a direction of vanishing exponents for its excluded constants that keeps every
cycle row of the catalytic graph (`_testhelper_cycle_basis`) with every other constant
fixed. A plain binding's equilibrium constant must fall (exponent at least 1, after
scaling); any other excluded edge's exponent is free.
"""
function _testhelper_lp_existence(am, c::Int)
    cm = EnzymeRates.catalytic(am)
    excluded = [k for k in eachindex(cm.edges) if _testhelper_excludes(am, k, c)]
    plain = [k for k in excluded if EnzymeRates._is_plain_binding(cm, cm.edges[k])]
    isempty(plain) && return true
    rows = _testhelper_cycle_basis((forms = cm.forms, steps = cm.edges))
    isempty(rows) && return true
    Z = Rational{BigInt}.(permutedims(reduce(hcat, rows)))
    AP, AF = Z[:, plain], Z[:, setdiff(excluded, plain)]
    _testhelper_feasible([AP AF -AF], -AP * ones(Rational{BigInt}, length(plain)))
end

"""
Whether the test's own certified direction for conformation `c` of `am` exists and keeps
every cycle row: the component potentials of `_testhelper_exclusion` (free E's at 0, every
other source at 1, arcs rising by 1) give each excluded edge the exponent of its bound
side over its unbound side; each arc's must be at least 1 and every row of
`_testhelper_cycle_basis` must sum to 0.
"""
function _testhelper_certified(am, c::Int)
    ex = _testhelper_exclusion(am, c)
    ex.failure === :none || return false
    cm = EnzymeRates.catalytic(am)
    start = [x == ex.component[1] ? big(0.0) : big(1.0) for x in 1:maximum(ex.component)]
    φ = _testhelper_component_potentials(ex, start, ones(BigFloat, length(ex.arcs)))
    a = [_testhelper_excludes(am, k, c) ? φ[ex.component[e.to]] - φ[ex.component[e.from]] :
         big(0.0) for (k, e) in enumerate(cm.edges)]
    rows = _testhelper_cycle_basis((forms = cm.forms, steps = cm.edges))
    all(arc -> a[arc[3]] >= 1, ex.arcs) && all(z -> iszero(sum(z .* a)), rows)
end

"""The uni-uni grid of admission's existence test: `E + S ⇌ E(S)` or `<-->`,
`E(S) <--> E(P)` and `E + P ⇌ E(P)` or `<-->`, each step with each of `:OnlyA`,
`:EqualAI` and `:NonequalAI`, 4 × 27 = 108 allosteric mechanisms, built unguarded."""
function _testhelper_uniuni_grid()
    tags = (:OnlyA, :EqualAI, :NonequalAI)
    [EnzymeRates.AllostericMechanism(EnzymeRates._raw_mechanism(Meta.parse("""
         begin
             substrates: S
             products: P
             steps: begin
                 E + S $s E(S) :: $tS
                 E(S) <--> E(P) :: $tC
                 E + P $p E(P) :: $tP
             end
         end"""), true))
     for s in ("⇌", "<-->") for p in ("⇌", "<-->")
     for tS in tags for tC in tags for tP in tags]
end

"""
The concerted network of `n` catalytic subunits: each subunit moves through its
conformation's law network (`_testhelper_network`) on its own, and the oligomer flips only
with every subunit at free E, at the single-subunit flip ratio. Its nodes are
`(state, conformation)`, a state numbering one assignment of the law network's nodes to
the subunits.
"""
function _testhelper_protomer_network(m, point, concs::AbstractDict; n::Int,
                                      speeds::Matrix{BigFloat})
    single = _testhelper_network(m, point, concs; speeds)
    conformation(u) = single.nodes[u][2]
    states = Vector{Int}[]
    for c in sort!(unique(last.(single.nodes)))
        members = [u for u in eachindex(single.nodes) if conformation(u) == c]
        append!(states, [collect(t) for t in Iterators.product(ntuple(_ -> members, n)...)])
    end
    index = Dict(s => q for (q, s) in enumerate(states))
    links = _testhelper_Link[]
    for (q, s) in enumerate(states), p in 1:n, l in single.links
        conformation(l.i) == conformation(l.j) && l.i == s[p] || continue
        t = copy(s)
        t[p] = l.j
        push!(links, _testhelper_Link(q, index[t], l.f, l.r, l.w, l.fast))
    end
    for l in single.links
        conformation(l.i) == conformation(l.j) && continue
        push!(links, _testhelper_Link(index[fill(l.i, n)], index[fill(l.j, n)], l.f, l.r,
                                      0, true))
    end
    _testhelper_Network([(q, conformation(s[1])) for (q, s) in enumerate(states)], links)
end

"""The steady-state turnover per active site of the concerted `n`-subunit network of `m`
(`_testhelper_protomer_network`), its rapid-equilibrium edges at `Λ` times `speeds`."""
function _testhelper_protomer_rate(m, point, concs::AbstractDict; n::Int,
                                   speeds::Matrix{BigFloat} = fill(big(1.0), size(point.K)),
                                   Λ::BigFloat = big(10.0)^30)
    network = _testhelper_protomer_network(m, point, concs; n, speeds)
    _testhelper_steady_rate(network; Λ).v / n
end

"""The guard's class-3 verdict (`_check_face`) on each face of `m`, in `_faces` order."""
function _testhelper_face_verdicts(m)
    mech = _testhelper_lift(m)
    net = EnzymeRates._law_network(mech)
    faces = EnzymeRates._faces(EnzymeRates.catalytic(mech))
    [EnzymeRates._check_face(mech, net, faces, i) for i in eachindex(faces)]
end

"""The guard's class-3 fault on each face of `m` (`_testhelper_face_verdicts`)."""
_testhelper_face_faults(m) = [v.fault for v in _testhelper_face_verdicts(m)]

"""Whether the guard's fault on a face agrees with the RE-speed oracle's verdict there: an
accepted face exactly when the oracle agrees, a rejected one exactly when it reads a speed
pattern, a second closed class, divergence, no law value or a mismatch."""
_testhelper_agrees(fault::EnzymeRates._Fault, verdict::Symbol) =
    fault === EnzymeRates._NO_FAULT ? verdict === :agree :
        verdict in (:pattern, :nonunique, :diverges, :law_nan, :mismatch)

"""Whether the guard's message for the class-3 verdict `v` on `m` names the face and, when
`v` leads with an edge, that edge's step line."""
function _testhelper_names_its_step(m, v::EnzymeRates._Verdict)
    mech = _testhelper_lift(m)
    msg = EnzymeRates._guard_message(mech, v)
    line = v.edge == 0 ? "" : EnzymeRates._step_text(EnzymeRates.catalytic(mech), v.edge)
    occursin(": at zero ", msg) && occursin(line, msg)
end

"""The faults class 3 reports on a face (`_check_face`): the faces the RE-speed oracle
judges."""
const _testhelper_class3_faults = (EnzymeRates._NONUNIQUE, EnzymeRates._BOTTOMLESS,
    EnzymeRates._TRAP, EnzymeRates._SPLIT, EnzymeRates._LAW_MULTI, EnzymeRates._MULTILEVEL)

"""
The guard's fault on each face of the allosteric mechanism `am`, in the order of `_faces`
(and `_testhelper_faces`), after its interior checks pass: class 3's where `_check_face`
rejects the face, else that of admission's face condition (`_admission_face`).
"""
function _testhelper_admission_faults(am)
    eps, law = EnzymeRates._eps_network(am), EnzymeRates._law_network(am)
    faces = EnzymeRates._faces(am.catalytic)
    map(eachindex(faces)) do i
        v = EnzymeRates._check_face(am, law, faces, i)
        v.fault === EnzymeRates._NO_FAULT ?
            EnzymeRates._admission_face(am, eps, law, faces, i).fault : v.fault
    end
end

"""Whether admission's `fault` on a face agrees with the ε-path oracle's `verdict` there:
an accepted face converges, a rejected one does not."""
_testhelper_admission_agrees(fault::EnzymeRates._Fault, verdict::Symbol) =
    fault === EnzymeRates._NO_FAULT ? verdict === :converges :
    verdict in (:no_path, :del_nonunique, :nonunique, :direction_dependent, :wrong_limit)

"""
Whether the guard's fault on each face of the allosteric mechanism `am`
(`_testhelper_admission_faults`) agrees with the oracle that judges it: a class-3 fault
with the RE-speed oracle (`_testhelper_agrees`), an admission fault with the ε-path oracle
(`_testhelper_admission_agrees`), and an accepted face with both. The RE-speed oracle
compares with its own deletion limit (`law = :limit`): the package cannot derive the law
of a mechanism the guard rejects. One `Bool` per face.
"""
function _testhelper_admission_agreement(am; rng)
    faults = _testhelper_admission_faults(am)
    faces = _testhelper_faces(am)
    respeed = _testhelper_respeed_verdicts(am; rng, law = :limit)
    map(eachindex(faults)) do i
        f = faults[i]
        f in _testhelper_class3_faults && return _testhelper_agrees(f, respeed[i])
        epath = _testhelper_admission_agrees(f, _testhelper_epath_verdict(am; rng,
                                                                          face = faces[i]))
        f === EnzymeRates._NO_FAULT ? epath && respeed[i] === :agree : epath
    end
end

"""
Check that the rate law of `m`, with one catalytic subunit, equals to 1e-20 the exact
rapid-equilibrium limit of the oracle's own law network (`_testhelper_network`,
`_testhelper_re_limit`), which deletes each conformation's excluded steps and prunes from
free E by itself, at `n_points` consistent points (`_testhelper_consistent_point`). A
compiled `m` is read through its derivation entry, which runs the guard; a concrete `m` is
read without it (`_testhelper_law_value`).
"""
function _testhelper_check_law_against_limit(m; n_points::Int = 3, seed::Int = 1)
    sub = _testhelper_single_subunit(m)
    am = _testhelper_lift(sub)
    fitted = EnzymeRates._dependent_param_exprs(
        sub isa EnzymeRates.AbstractEnzymeMechanism ? typeof(sub) : sub)[2]
    rng = Random.MersenneTwister(seed)
    setprecision(BigFloat, 256) do
        for _ in 1:n_points
            point, concs, params = _testhelper_law_inputs(am, fitted, rng)
            speeds = _testhelper_draw_speeds(am, rng)
            exact = _testhelper_re_limit(_testhelper_network(am, point, concs; speeds))
            @test abs(_testhelper_law_value(sub, params, concs) - exact) <=
                  big(10)^-20 * abs(exact)
        end
    end
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
        analytical_peak_turnover_fn=p -> p.k3f * p.k4f / (p.k3f + p.k4f),
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
        analytical_peak_turnover_fn=p -> p.k3f,
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
        analytical_peak_turnover_fn=p -> p.k2f * p.k4f / (p.k2f + p.k4f),
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
            E + S1 ⇌ E(S1)
            E + S2 ⇌ E(S2)
            E(S1) + S2 ⇌ E(S1, S2)
            E(S2) + S1 <--> E(S1, S2)
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
        analytical_peak_turnover_fn=p -> p.k2f,
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
        analytical_peak_turnover_fn=p -> p.k2f,
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
        analytical_peak_turnover_fn=p -> p.k2f,
        expected_factored_num=
        "k_ERS_to_EPR * R * S / (K_ERS_to_ER_S * K_ER_to_E_R)" *
        " - k_EPR_to_ERS * P * R / (K_EPR_to_ER_P * K_ER_to_E_R)",
        expected_factored_denom=
        "1 + R / K_ER_to_E_R + P * R / (K_EPR_to_ER_P * K_ER_to_E_R)" *
        " + R * S / (K_ERS_to_ER_S * K_ER_to_E_R)",
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
        analytical_peak_turnover_fn=p -> max(p.k2f, p.k5f),
        expected_factored_num=
        "k_ES_to_EP * S / K_ES_to_E_S" *
        " + k_ERS_to_EPR * R * S / (K_ERS_to_ER_S * K_ER_to_E_R)" *
        " - (k_EP_to_ES * P / K_EP_to_E_P" *
        " + k_EPR_to_ERS * P * R / (K_EP_to_E_P * K_ER_to_E_R))",
        expected_factored_denom=
        "1 + P / K_EP_to_E_P + R / K_ER_to_E_R + S / K_ES_to_E_S" *
        " + P * R / (K_EP_to_E_P * K_ER_to_E_R) + R * S / (K_ERS_to_ER_S * K_ER_to_E_R)",
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
        "k_ES_to_EP * S / K_ES_to_E_S" *
        " + k_EAS_to_EAP * A * S / (K_EAS_to_EA_S * K_EA_to_E_A)" *
        " - (k_EP_to_ES * P / K_EP_to_E_P" *
        " + k_EAP_to_EAS * A * P / (K_EAP_to_EA_P * K_EA_to_E_A))",
        expected_factored_denom=
        "1 + A / K_EA_to_E_A + I / K_EI_to_E_I + P / K_EP_to_E_P + S / K_ES_to_E_S" *
        " + A * P / (K_EAP_to_EA_P * K_EA_to_E_A) + A * S / (K_EAS_to_EA_S * K_EA_to_E_A)",
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
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S" *
        " + I1 * I2 / (K_EI1_to_E_I1 * K_EI2_to_E_I2)" *
        " + I1 * P / (K_EI1_to_E_I1 * K_EP_to_E_P)" *
        " + I1 * S / (K_EI1_to_E_I1 * K_ES_to_E_S)" *
        " + I2 * P / (K_EI2_to_E_I2 * K_EP_to_E_P)" *
        " + I2 * S / (K_EI2_to_E_I2 * K_ES_to_E_S)" *
        " + I1 * I2 * P / (K_EI1_to_E_I1 * K_EI2_to_E_I2 * K_EP_to_E_P)" *
        " + I1 * I2 * S / (K_EI1_to_E_I1 * K_EI2_to_E_I2 * K_ES_to_E_S)",
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
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S" *
        " + I1 * P / (K_EI1_to_E_I1 * K_EP_to_E_P)" *
        " + I1 * S / (K_EI1_to_E_I1 * K_ES_to_E_S)",
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
        "1 + I1 / K_EI1_to_E_I1 + I2 / K_EI2_to_E_I2 + P / K_EP_to_E_P + S / K_ES_to_E_S" *
        " + I1 * P / (K_EI1_to_E_I1 * K_EP_to_E_P)" *
        " + I1 * S / (K_EI1_to_E_I1 * K_ES_to_E_S)" *
        " + I2 * P / (K_EI2_to_E_I2 * K_EP_to_E_P)" *
        " + I2 * S / (K_EI2_to_E_I2 * K_ES_to_E_S)",
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
        analytical_peak_turnover_fn = p -> p.k5f,
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
        analytical_peak_turnover_fn = p -> p.k6f,
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
        analytical_peak_turnover_fn = p -> p.k5f,
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
        analytical_peak_turnover_fn=nothing,      # cat is :NonequalAI → kcat L-dependent
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
        analytical_peak_turnover_fn=_testhelper_m_OnlyA_prod_peak_turnover,
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
            (E + NAD <--> E(NAD),
             E(Lactate) + NAD <--> E(Lactate, NAD),
             E(Pyruvate) + NAD <--> E(NAD, Pyruvate)) :: EqualAI
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
        name="LDH i-state NonequalAI 6-group, steady-state NAD group",
        mechanism=m,
        raw=raw,
        metabolite_names=[:NADH, :Pyruvate, :Lactate, :NAD],
        expected_n_states=9, expected_n_steps=13,
        expected_n_haldane_constraints=1,
        expected_n_wegscheider_constraints=2, expected_n_independent_params=10,
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
