# ABOUTME: Plain-data mechanism populations (module Pop) shared by the derivation probes.
# ABOUTME: Mirrors EnzymeRates mechanisms plus the package's outputs; stdlib only.
module Pop

using Serialization

const DATA_DIR = joinpath(@__DIR__, "data")

"""A metabolite as a step or form holds it. `role` is `:S` substrate, `:P` product,
`:Inh` competitive inhibitor (a standalone inhibitor `I` or a reactant's copy `A::Inh`,
which keeps the reactant's name and shares its concentration) or `:Reg` allosteric
regulator."""
struct PMet
    name::Symbol
    role::Symbol
end

"""An enzyme form. `bound` is sorted as the package sorts it (by name, an inhibitor copy
after the reactant of the same name). `name` is the package's rendered name."""
struct PSpecies
    bound::Vector{PMet}
    conformation::Symbol
    res_added::Vector{Symbol}
    res_subtracted::Vector{Symbol}
    name::Symbol
end

"""One step in the package's STORED direction: `from` takes up `consumed` and gives off
`released` to become `to` (indices into `PMech.species`). `re` marks rapid equilibrium."""
struct PStep
    from::Int
    to::Int
    consumed::Vector{PMet}
    released::Vector{PMet}
    re::Bool
end

"""An allosteric regulatory site: ligands with their per-ligand state tags."""
struct PSite
    ligands::Vector{Symbol}
    multiplicity::Int
    states::Vector{Symbol}
end

"""A monomial as `symbol => exponent` pairs sorted by symbol; a polynomial as
`monomial => coefficient` pairs, sorted by the monomial's text."""
const PMono = Vector{Pair{Symbol, Int}}
const PPoly = Vector{Pair{PMono, Rational{Int}}}

"""One numeric sample of the package. `params` are the Reduced parameters passed to
`rate_equation` (fitted..., Keq, E_total = 1); `dep_values` every dependent the package
derives from them (its own expressions, evaluated in Float64); `v` the Reduced rate.
`full_params`/`v_full` (plain mechanisms only) are random raw group constants
(`parameters(em, Full)`, E_total = 1) and the Full-mode rate at the same `concs`. `kcat` is
`_kcat_forward(em, params)`; NaN when not computed or when it threw."""
struct PSample
    concs::Vector{Pair{Symbol, Float64}}
    params::Vector{Pair{Symbol, Float64}}
    dep_values::Vector{Pair{Symbol, Float64}}
    v::Float64
    full_params::Vector{Pair{Symbol, Float64}}
    v_full::Float64
    kcat::Float64
end

"""One mechanism with its structure and the package's outputs. Steps are flat, in the
package's stored order (group by group); `groups[g]` lists the flat indices of kinetic
group `g`. See README.md for every field."""
Base.@kwdef struct PMech
    id::String
    source::String
    family::Vector{Symbol} = Symbol[]
    substrates::Vector{Symbol}
    products::Vector{Symbol}
    inhibitors::Vector{Symbol} = Symbol[]
    regulators::Vector{Symbol} = Symbol[]
    metabolites::Vector{Symbol} = Symbol[]
    species::Vector{PSpecies}
    steps::Vector{PStep}
    groups::Vector{Vector{Int}}
    allosteric::Bool = false
    cat_multiplicity::Int = 1
    group_tags::Vector{Symbol} = Symbol[]
    sites::Vector{PSite} = PSite[]
    ok::Bool = true
    error::String = ""
    stage_errors::Vector{Pair{Symbol, String}} = Pair{Symbol, String}[]
    step_names::Vector{Vector{Symbol}} = Vector{Symbol}[]
    step_names_A::Vector{Vector{Symbol}} = Vector{Symbol}[]
    step_names_I::Vector{Vector{Symbol}} = Vector{Symbol}[]
    step_names_I_pkg::Vector{Vector{Symbol}} = Vector{Symbol}[]
    inactive_kept::Vector{Bool} = Bool[]
    inactive_reversed::Vector{Bool} = Bool[]
    fitted_params::Vector{Symbol} = Symbol[]
    full_params::Vector{Symbol} = Symbol[]
    dependents::Vector{Pair{Symbol, String}} = Pair{Symbol, String}[]
    dependent_monomials::Vector{Pair{Symbol, Vector{Pair{Symbol, Rational{Int}}}}} =
        Pair{Symbol, Vector{Pair{Symbol, Rational{Int}}}}[]
    rate_equation_string::String = ""
    eq_hash::UInt64 = 0
    independent_param_count::Int = -1
    pass2_rename::Vector{Pair{Symbol, Symbol}} = Pair{Symbol, Symbol}[]
    pass2_rename_A::Vector{Pair{Symbol, Symbol}} = Pair{Symbol, Symbol}[]
    pass2_rename_I::Vector{Pair{Symbol, Symbol}} = Pair{Symbol, Symbol}[]
    pass2_changes_polys::Bool = false
    raw_num::PPoly = PPoly()
    raw_den::PPoly = PPoly()
    raw_dfree::PPoly = PPoly()
    pkg_num::PPoly = PPoly()
    pkg_den::PPoly = PPoly()
    pkg_dfree::PPoly = PPoly()
    state_polys_A::Vector{PPoly} = PPoly[]
    state_polys_I::Vector{PPoly} = PPoly[]
    samples::Vector{PSample} = PSample[]
    eq_complexity::Int = -1
    ctor_ns::Float64 = NaN
    first_call_s::Float64 = NaN
    eval_ns::Float64 = NaN
    eval_allocs::Int = -1
end

"""The population stored as `data/<name>.jls`."""
load_population(name::AbstractString) =
    deserialize(joinpath(DATA_DIR, string(name, ".jls")))::Vector{PMech}

"""Write `mechs` as population `name` (or to an explicit `path`)."""
save_population(name::AbstractString, mechs::Vector{PMech};
                path = joinpath(DATA_DIR, string(name, ".jls"))) =
    (serialize(path, mechs); path)

"""The populations present in `data/`."""
populations() = sort!([splitext(f)[1] for f in readdir(DATA_DIR)
                       if endswith(f, ".jls") && !occursin(".part", f)])

"""Each kinetic group of every flat step."""
group_of_step(m::PMech) = (g = zeros(Int, length(m.steps));
                           for (k, grp) in enumerate(m.groups), j in grp; g[j] = k; end; g)

"""Whether `s` is a binding in the package's sense: one metabolite consumed, none
released (plain or fused)."""
is_binding(s::PStep) = length(s.consumed) == 1 && isempty(s.released)

"""Whether the form holds no metabolite and no residual."""
is_free(sp::PSpecies) =
    isempty(sp.bound) && isempty(sp.res_added) && isempty(sp.res_subtracted)

"""Value of a dependent's monomial (`symbol => exponent`, `:Keq` included) at `vals`."""
eval_monomial(mono, vals) = prod((Float64(vals[s])^Float64(e) for (s, e) in mono); init = 1.0)

"""Number of independent RE cycles of a plain mechanism (it reads `step_names`) and how
many of them the kinetic grouping does not close identically: around a fundamental cycle
of the RE subgraph the log-weights must return to zero, and each step adds ∓log K of its
group's name (− for a binding's dissociation constant, + otherwise); a cycle whose name
coefficients cancel is closed by the grouping itself, any other carries a Wegscheider
relation that Full mode does not impose."""
function re_cycle_closure(m)
    n = length(m.species)
    adj = [Tuple{Int, Int, Int}[] for _ in 1:n]          # (neighbor, step, sign)
    for (j, s) in enumerate(m.steps)
        s.re || continue
        push!(adj[s.from], (s.to, j, 1)); push!(adj[s.to], (s.from, j, -1))
    end
    coef(j) = is_binding(m.steps[j]) ? -1 : 1
    pot = Vector{Union{Nothing, Dict{Symbol, Int}}}(nothing, n)
    tree = Set{Int}()
    for r in 1:n
        (pot[r] === nothing && !isempty(adj[r])) || continue
        pot[r] = Dict{Symbol, Int}(); queue = [r]
        while !isempty(queue)
            u = popfirst!(queue)
            for (v, j, sg) in adj[u]
                pot[v] === nothing || continue
                d = copy(pot[u]); k = m.step_names[j][1]
                d[k] = get(d, k, 0) + sg * coef(j)
                pot[v] = d; push!(tree, j); push!(queue, v)
            end
        end
    end
    cycles = 0; open = 0
    for (j, s) in enumerate(m.steps)
        (s.re && j ∉ tree) || continue
        cycles += 1
        d = mergewith(+, pot[s.to], Dict(k => -v for (k, v) in pot[s.from]))
        k = m.step_names[j][1]
        d[k] = get(d, k, 0) - coef(j)
        any(!iszero, values(d)) && (open += 1)
    end
    cycles, open
end

"""A sample's concentrations, Reduced parameters, dependents, or Full parameters as
a `Dict`."""
concs(s::PSample) = Dict(s.concs)
params(s::PSample) = Dict(s.params)
all_values(s::PSample) = merge(Dict(s.params), Dict(s.dep_values))
full_params(s::PSample) = Dict(s.full_params)

"""A step as the package's DSL would print it, in its stored direction."""
function step_text(m::PMech, s::PStep)
    side(f, ms) = join([String(m.species[f].name); [String(x.name) * (x.role === :Inh ?
                        "::Inh" : "") for x in ms]], " + ")
    side(s.from, s.consumed) * (s.re ? " ⇌ " : " <--> ") * side(s.to, s.released)
end

end # module
