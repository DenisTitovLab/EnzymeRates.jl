# ABOUTME: Haldane/Wegscheider thermodynamic constraints for enzyme mechanisms.
# ABOUTME: Finds cycles, selects dependent params, builds @generated rate-eq preambles.

"""
Haldane and Wegscheider thermodynamic constraints for enzyme mechanisms,
plus preamble building helpers for @generated rate equation bodies.

Identifies thermodynamic cycles in the mechanism graph (via null-space of the
enzyme incidence matrix), classifies them as Haldane (net reaction) or
Wegscheider (internal loop), and performs Gaussian elimination to express a
minimal set of dependent rate constants in terms of the independent ones and
the equilibrium constant Keq.
"""

# ─── Shared Helpers ──────────────────────────────────────────────

"""
Collect raw parameter symbols (K_i for RE, k_if/k_ir for SS) for the
representative step of each kinetic group, in step-order. Routes
Symbol production through the `name(p::Parameter, m)` chokepoint via
`_enumerate_parameters_full`.
"""
function _raw_param_symbols(m::Mechanism)
    Symbol[name(p, m) for p in _enumerate_parameters_full(m)]
end

_raw_param_symbols(m::EnzymeMechanism) = _raw_param_symbols(Mechanism(m))

"""
For each step in `m` (in flat-iteration order), yield the Parameter
instances that govern that step: `[Kd|Kiso]` for an RE step, `[Kon|Kfor,
Koff|Krev]` for an SS step. Each Parameter is anchored on the original
step (not the rep), so `name(p, m)` renders to the rep's structural
Symbol via the value-context chokepoint, collapsing kinetic-group members.
"""
function _step_parameters(m::Mechanism)
    out = Vector{Vector{Parameter}}()
    for (s, _) in _flat_steps(m)
        params = is_equilibrium(s) ?
            Parameter[is_binding(s) ? Kd(s, :None) : Kiso(s, :None)] :
            Parameter[is_binding(s) ? Kon(s, :None)  : Kfor(s, :None),
                      is_binding(s) ? Koff(s, :None) : Krev(s, :None)]
        push!(out, params)
    end
    out
end

# ─── Structural primacy: free-enzyme set + step priority ─────────

"""
Set of enzyme-form names that are NOT the RHS of any canonical RE
binding step `F + met… ⇌ F_bound`. Walks `Mechanism.steps` directly:
for each RE binding step, the canonical form puts the bound metabolite
on `to_species`, so `to_species`'s name is excluded from the free set.
Iso steps don't determine binding state. SS steps' direction is not
canonicalized so they don't participate.

A form that carries bound metabolites but has no binding-in step in this
graph is also excluded: an inactive-conformation graph
(`_state_mechanism(am, :I)`) drops each `:OnlyA` binding step, so the ligand's
downstream complex is reached only by conformational flip or reverse catalysis
and loses its incoming binding edge — yet it is still a bound form, never free
enzyme. Excluding it keeps a kinetic group's naming representative invariant
across the active/inactive split (an `:EqualAI` group renders one shared
Symbol), which is a no-op on canonical mechanisms where every bound form has a
binding-in step.

Shared by the kinetic-group name representative and the Haldane
elimination pivot.
"""
function _free_enz_set(m::Union{Mechanism, AllostericMechanism})
    enz_names = Set{Symbol}()
    for group in steps(m), s in group
        push!(enz_names, name(from_species(s)))
        push!(enz_names, name(to_species(s)))
    end
    free_enz_set = copy(enz_names)
    for group in steps(m), s in group
        is_equilibrium(s) || continue
        is_binding(s) || continue
        # Canonical: bound metabolite resides on to_species. The from-side
        # is the "free + met" reactant; the to-side is the bound form.
        delete!(free_enz_set, name(to_species(s)))
    end
    bound_in = Set{Symbol}(name(to_species(s))
                           for group in steps(m) for s in group if is_binding(s))
    for group in steps(m), s in group
        for sp in (from_species(s), to_species(s))
            isempty(bound(sp)) || name(sp) in bound_in ||
                delete!(free_enz_set, name(sp))
        end
    end
    free_enz_set
end

"""
Structural primacy base score for a step (lower = more primary / less
eliminable). Free-enzyme RE binding (-1) < free-enzyme SS binding (0) <
non-free metabolite step (10) < internal isomerization (20). Shared by the
kinetic-group name representative (argmin) and the Haldane elimination pivot
(argmax, which adds a +0/+1 forward/reverse offset per rate constant).
"""
function _step_priority(s::Step, free_enz_set::Set{Symbol})
    has_met = is_binding(s)
    is_free = (name(from_species(s)) in free_enz_set) ||
              (name(to_species(s))   in free_enz_set)
    is_equilibrium(s) && has_met && is_free && return -1
    return !has_met ? 20 : is_free ? 0 : 10
end

"""
Total lexical tiebreak for two distinct steps in the same kinetic group:
species pair + bound metabolite + RE/SS flag.
"""
_step_lex_key(s::Step) =
    (String(name(from_species(s))), String(name(to_species(s))),
     String(bound_metabolite(s) === nothing ? "" : name(bound_metabolite(s))),
     is_equilibrium(s))

"""
Kinetic-group naming representative: the structurally-primary step
(`argmin _step_priority`), with a deterministic lexical tiebreak.
"""
_group_rep(group::Vector{Step}, free_enz_set::Set{Symbol}) =
    argmin(s -> (_step_priority(s, free_enz_set), _step_lex_key(s)), group)

# ─── Thermodynamic Constraint Infrastructure ─────────────────────

# Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
# free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
# pivot at reduced-matrix row i) plus the reduced matrix R. Used by
# `_integer_nullspace` (nullspace basis).
function _rref_partition(A::Matrix{Int})
    m, n = size(A)
    R = Matrix{Rational{BigInt}}(A)
    pivot_cols = Int[]
    row = 1
    for col in 1:n
        piv = findfirst(r -> R[r, col] != 0, row:m)
        piv === nothing && continue
        piv += row - 1
        R[row, :], R[piv, :] = R[piv, :], R[row, :]
        R[row, :] ./= R[row, col]
        for r in 1:m
            r != row && R[r, col] != 0 && (R[r, :] .-= R[r, col] .* R[row, :])
        end
        push!(pivot_cols, col); row += 1
    end
    free_cols = setdiff(1:n, pivot_cols)
    return pivot_cols, free_cols, R
end

function _integer_nullspace(A::Matrix{Int})
    n = size(A, 2)
    pivot_cols, free_cols, R = _rref_partition(A)
    isempty(free_cols) && return zeros(Int, n, 0)
    NS = zeros(Rational{BigInt}, n, length(free_cols))
    for (k, fc) in enumerate(free_cols)
        NS[fc, k] = 1
        for (r, pc) in enumerate(pivot_cols); NS[pc, k] = -R[r, fc]; end
    end
    result = zeros(Int, n, length(free_cols))
    for k in axes(result, 2)
        col = @view NS[:, k]
        l = lcm(denominator.(col)...)
        result[:, k] .= Int.(col .* l)
        g = gcd(abs.(result[:, k])...)
        g > 0 && (result[:, k] .÷= g)
        fnz = findfirst(!=(0), result[:, k])
        fnz !== nothing && result[fnz, k] < 0 && (result[:, k] .*= -1)
    end
    result
end

function _thermodynamic_constraints(mech::Mechanism)
    flat = _flat_steps(mech)
    enz_names = collect(_enumerate_species_names(mech))
    enz_name_to_idx = Dict(n => i for (i, n) in enumerate(enz_names))
    nsteps = length(flat)
    met_names = Symbol[name(metabolite(ra)) for ra in reactants(mech.reaction)]
    subs_species = Symbol[name(s) for s in substrates(mech.reaction)]
    prods_species = Symbol[name(p) for p in products(mech.reaction)]

    # Enzyme incidence matrix
    B = zeros(Int, length(enz_names), nsteps)
    for (j, (s, _)) in enumerate(flat)
        i_from = enz_name_to_idx[name(from_species(s))]
        i_to   = enz_name_to_idx[name(to_species(s))]
        B[i_from, j] -= 1
        B[i_to,   j] += 1
    end

    # Stoichiometry matrix (rows = metabolites, cols = steps). A
    # metabolite gets its stoichiometry solely from the canonical
    # reaction tuple via `_step_sides(s)`: m_lhs contributes -1 (consumed
    # from the free pool), m_rhs contributes +1 (produced).
    #
    # Iso steps carry no free-pool metabolite — their bound content is
    # encoded in the enzyme-form identity — so `_step_sides` returns empty
    # metabolite lists and they contribute zero. Do NOT add a
    # from_bound/to_bound diff for iso steps: that double-counts
    # metabolites already accounted for by the binding/release steps and
    # inflates the cycle's net change (e.g. 1/Keq -> 1/Keq^2).
    met_idx = Dict(n => i for (i, n) in enumerate(met_names))
    stoich_mat = zeros(Int, length(met_names), nsteps)
    for (j, (s, _)) in enumerate(flat)
        _, _, m_lhs, m_rhs = _step_sides(s)
        for m in m_lhs
            haskey(met_idx, m) && (stoich_mat[met_idx[m], j] -= 1)
        end
        for m in m_rhs
            haskey(met_idx, m) && (stoich_mat[met_idx[m], j] += 1)
        end
    end

    NS = _integer_nullspace(B)
    nc = size(NS, 2)
    nc == 0 && return zeros(Int, 0, size(B, 2)), Int[]

    nu_net = zeros(Int, length(met_names))
    for nm in subs_species
        nu_net[met_idx[nm]] -= 1
    end
    for nm in prods_species
        nu_net[met_idx[nm]] += 1
    end

    # Classify each null-space cycle as Haldane (proportional to the
    # net reaction → contributes log(Keq)) or Wegscheider (closed
    # cycle, zero net change). Errors on cycles that touch metabolites
    # but aren't proportional to the net reaction.
    function classify_cycle(nu_cycle, i)
        all(nu_cycle .== 0) && return 0
        c = nothing
        for j in eachindex(nu_cycle)
            if nu_net[j] == 0
                nu_cycle[j] != 0 && error(
                    "Cycle $i produces metabolite " *
                    "change not proportional to " *
                    "net reaction"
                )
            else
                c_j = nu_cycle[j] // nu_net[j]
                if c === nothing
                    c = c_j
                elseif c_j != c
                    error(
                        "Cycle $i produces metabolite " *
                        "change not proportional to " *
                        "net reaction"
                    )
                end
            end
        end
        err = "Cycle $i produces metabolite change " *
              "not proportional to net reaction"
        c !== nothing && denominator(c) == 1 ? Int(c) : error(err)
    end

    C = NS'
    rhs_coeffs = [classify_cycle(stoich_mat * C[i, :], i) for i in 1:nc]
    return C, rhs_coeffs
end

# Walk Mechanism.steps; emit distinct enzyme-form Symbol names in
# step-walk order. Used by _thermodynamic_constraints and friends.
function _enumerate_species_names(mech::Mechanism)
    seen = Symbol[]
    for group in steps(mech), s in group
        for sp in (from_species(s), to_species(s))
            nm = name(sp)
            nm in seen || push!(seen, nm)
        end
    end
    seen
end

"""
    _dependent_param_exprs(M::Type{<:EnzymeMechanism}) → (dep_exprs, indep_params)

Select dependent parameters and build substitution expressions for the
Haldane / Wegscheider thermodynamic constraints. Steps in the same
kinetic group share parameters: their cycle-incidence columns are merged
into the representative step's column before Gaussian elimination, so
`dep_exprs` and `indep_params` are keyed only on representatives.

The wrapper calls `_build_wegscheider_rename_map(M)` to obtain the
rename map for absorbed single-symbol Wegscheider RE ties and forwards
to `_dependent_param_exprs_kernel`.
"""
function _dependent_param_exprs(M::Type{<:EnzymeMechanism})
    rename = _build_wegscheider_rename_map(M)
    dep_exprs, indep = _dependent_param_exprs_kernel(M, rename)
    # Filter Pass-2-absorbed symbols out of indep. Pass 2 of
    # `_build_wegscheider_rename_map` adds entries like `K_P_E => K_S_E`
    # when a Wegscheider tie collapses two binding-K group reps to the
    # same name. After the merge, the absorbed symbol doesn't appear in
    # the v polynomial — its column has been folded into the target.
    # But the absorbed symbol is still a kinetic-group rep in the
    # mechanism, so `_raw_param_symbols` emits it and the kernel keeps it
    # in `indep`. Without this filter, `fitted_params` exposes a fittable
    # dummy dimension that doesn't affect the loss, and finite-restart
    # convergence suffers (the same rate equation can land at noticeably
    # different fitted losses depending on which absorbed symbol got
    # the dummy slot).
    indep = Tuple(p for p in indep if get(rename, p, p) == p)
    # Sort by name so `fitted_params` / the params destructuring is
    # content-canonical: two mechanisms with the same independent set (e.g.
    # graph-distinct but rate-equivalent ones) produce the identical rate
    # equation string and therefore the same dedup key.
    indep = Tuple(sort(collect(indep); by = string))
    return dep_exprs, indep
end

"""
Assemble the rational thermodynamic-constraint system for `mech`. Returns
`(A, rhs, columns, priority)`: `A` is the constraint matrix (rows = independent
Wegscheider/Haldane cycles, columns = parameters in `all_params` order), `rhs`
the per-row `log(Keq)` exponent, `columns` the ordered parameter symbols, and
`priority` the per-column pivot preference as an `(is_i_state, type)` tuple —
lexicographic, so an I-state column (`is_i_state = true`, set via the
`is_i_state` kwarg) always outranks an A-state/non-allosteric one, and within a
state `type` orders internal isomerizations > metabolite steps > free-enzyme
binding (higher scores are eliminated first, i.e. become dependent).
`_solve_dependent_set` consumes this. Split out from the kernel so the
allosteric derivation can stack per-state systems and reuse one solver.

Binding K's are Kd in the polynomial while cycle products use 1/Kd, so binding-K
column entries carry a sign flip on top of the cycle incidence. Non-representative
steps fold into their representative through the `name(p, mech)` chokepoint (plus
any Pass-2 single-symbol Wegscheider tie in `rename`) — equivalent to a
kinetic-group equality constraint.
"""
function _assemble_constraints(
    mech::Mechanism,
    rename::AbstractDict{Symbol, Symbol};
    step_params = _step_parameters(mech),
    all_params = _raw_param_symbols(mech),
    is_i_state::Bool = false,
)
    flat = _flat_steps(mech)
    free_enz_set = _free_enz_set(mech)

    C, rhs_coeffs = _thermodynamic_constraints(mech)
    nc = size(C, 1)
    nsteps = size(C, 2)

    columns = collect(all_params)
    sym_col = Dict(p => i for (i, p) in enumerate(columns))
    n_vars = length(columns)

    step_name(p::Parameter) = get(rename, name(p, mech), name(p, mech))

    binding_K_set = Set{Symbol}()
    for (j, (s, _)) in enumerate(flat)
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_K_set, step_name(step_params[j][1]))
    end

    A = zeros(Rational{BigInt}, nc, n_vars)
    rhs = Rational{BigInt}.(rhs_coeffs)
    for i in 1:nc, j in 1:nsteps
        C[i, j] == 0 && continue
        if is_equilibrium(flat[j][1])
            sym = step_name(step_params[j][1])
            sign_factor = sym in binding_K_set ? -1 : 1
            A[i, sym_col[sym]] += sign_factor * C[i, j]
        else
            kf = step_name(step_params[j][1])
            kr = step_name(step_params[j][2])
            A[i, sym_col[kf]] += C[i, j]
            A[i, sym_col[kr]] -= C[i, j]
        end
    end

    # Pivot priority: (is_I_state, type_priority). Lexicographic — an I-state column
    # outranks any A-state / non-allosteric column, so a cross-state affinity split
    # collapses onto the free A-side; within a state the `_step_priority` order holds.
    # No value is a never-pivot sentinel (that lives only in `_solve_dependent_set`).
    priority = fill((is_i_state, 0), n_vars)
    for j in 1:nsteps
        step = step_params[j][1].step
        base = _step_priority(step, free_enz_set)
        if is_equilibrium(flat[j][1])
            s = step_name(step_params[j][1])
            haskey(sym_col, s) && (priority[sym_col[s]] = (is_i_state, base))
        else
            for (offset, p) in enumerate(step_params[j])
                s = step_name(p)
                haskey(sym_col, s) &&
                    (priority[sym_col[s]] = (is_i_state, base + offset - 1))
            end
        end
    end
    return A, rhs, columns, priority
end

"""
    _onlya_haldane_violation(rxn, cat_steps, cat_allo_states)
        → Union{Nothing, String}

Return `nothing` when every catalytic thermodynamic cycle's Haldane relation
stays satisfiable under the `K_I = K_A/ε`, `ε → 0⁺` limit an `:OnlyA` binding
asserts; otherwise return a message naming the offending `:OnlyA` bindings.

A cycle's Haldane carries `∏ε_p/∏ε_s` on the inactive side. The `ε` are
independent, so that monomial can be held finite exactly when its exponents
carry both signs. All-same-sign drives it to `0` or `∞`, and only `k_I = 0` —
an `:OnlyA` catalytic tag — absorbs that.

The check graph drops `:OnlyA` chemical groups, so a cycle running through one
never appears and never reports a violation: that is the `k_I = 0` escape.
Bindings completing no cycle (competitive inhibitors, dead ends, regulator
sites) never enter a row and take no part. Both catalytic (Haldane) and
binding-only (Wegscheider, `rhs = 0`) cycle rows are inspected, so a one-sided
`:OnlyA` binding on a pure random-order binding square is caught.

The per-row sign test is a sufficient rejection condition, not a complete one:
it flags a cycle only when that single row's `:OnlyA` exponents are all one
sign. A genuinely-complete test asks whether the coupled `ε`-exponent system has
a strictly-positive nullspace vector (Stiemke feasibility). The two agree for
every random-order mechanism up to bi-bi; from ter-substrate up, a multi-cycle
coupled inconsistency can pass this per-row check. Such a mechanism is still
derived correctly — `:OnlyA` deletes the offending edges, so the rate law is
that of a consistent subgraph — so the gap is a checker-completeness contract
issue, not a wrong-equation one.

An RE binding carries the cycle exponent on its `Kd` column, already sign-flipped
against the cycle's `1/Kd` product, while an SS binding carries it on `Kon`
unflipped. Both encode the same `-C·log(Kd)`, so each column is normalized back
to the `ε` exponent before its sign is read; otherwise a cycle mixing the two
step kinds reads as same-sign and a balanced pair is rejected.

Builds a plain `Mechanism`; it must not call `_state_allo_mechanism`, which
would construct an `AllostericMechanism` and recurse.
"""
function _onlya_haldane_violation(rxn::EnzymeReaction,
                                  cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    keep = [g for g in eachindex(cat_steps)
            if !(cat_allo_states[g] === :OnlyA && is_iso(cat_steps[g][1]))]
    isempty(keep) && return nothing
    onlyA_steps = Set{Step}()
    for g in eachindex(cat_steps)
        cat_allo_states[g] === :OnlyA && is_binding(cat_steps[g][1]) &&
            union!(onlyA_steps, cat_steps[g])
    end
    isempty(onlyA_steps) && return nothing
    cm = Mechanism(rxn, [copy(cat_steps[g]) for g in keep])
    sp = _step_parameters(cm)
    A, _, columns, _ = _assemble_constraints(cm, Dict{Symbol, Symbol}();
                                             step_params = sp)
    sym_col = Dict(c => i for (i, c) in enumerate(columns))
    # Column → multiplier normalizing its entry to the ε exponent. Mirrors the
    # `sym in binding_K_set` sign rule of `_assemble_constraints` (line 366):
    # only an RE binding lands in that set, so only it is sign-flipped. The
    # multiplier is well defined per column because RE and SS bindings render to
    # different symbols (`K_S_E` vs `kon_S_E`) and so never share a column.
    onlyA_cols = Dict{Int, Int}()
    for (j, (s, _)) in enumerate(_flat_steps(cm))
        s in onlyA_steps || continue
        sym = name(sp[j][1], cm)
        haskey(sym_col, sym) &&
            (onlyA_cols[sym_col[sym]] = is_equilibrium(s) ? -1 : 1)
    end
    isempty(onlyA_cols) && return nothing
    for i in axes(A, 1)
        signs = Set{Int}()
        for (c, mult) in onlyA_cols
            A[i, c] == 0 || push!(signs, mult * A[i, c] > 0 ? 1 : -1)
        end
        isempty(signs) && continue
        length(signs) == 1 || continue
        offenders = sort!([string(columns[c]) for c in keys(onlyA_cols)
                           if A[i, c] != 0])
        return "an :OnlyA binding ($(join(offenders, ", "))) leaves a " *
               "thermodynamic (Haldane/Wegscheider) cycle unsatisfiable: the " *
               "inactive conformation cannot close that cycle at finite nonzero " *
               "affinity. Tag the cycle's chemical step :OnlyA, or tag an " *
               "opposing binding :OnlyA so the affinities diverge together."
    end
    nothing
end

"""
Solve an assembled constraint system `(A, rhs, columns, priority)` for a
dependent/independent parameter partition. Gaussian elimination with priority
pivoting picks, per constraint row, the highest-priority still-unused column as
the pivot (that column becomes dependent, expressed via the remaining columns);
every non-pivot column is independent (fitted). A row with no eligible pivot is
either redundant (`0 = 0`) or a thermodynamic contradiction (`0 = c·log Keq`),
which errors. Returns `(dep_exprs, indep)`. An empty system (no rows) yields no
dependents and all columns independent.
"""
function _solve_dependent_set(
    A::AbstractMatrix{Rational{BigInt}},
    rhs::AbstractVector{Rational{BigInt}},
    columns::AbstractVector{Symbol},
    priority::AbstractVector{Tuple{Bool, Int}},
)
    nc = size(A, 1)
    n_vars = length(columns)

    pivot_entries = Tuple{Int, Int}[]
    pivot_col_set = Set{Int}()
    wA, wrhs = copy(A), copy(rhs)
    for i in 1:nc
        best_col, best_pri = 0, (false, typemin(Int))
        for c in 1:n_vars
            c in pivot_col_set && continue
            wA[i, c] == 0 && continue
            if priority[c] > best_pri
                best_pri = priority[c]
                best_col = c
            end
        end
        if best_col == 0
            wrhs[i] == 0 && continue  # redundant constraint (0 = 0)
            error(
                "Thermodynamically contradictory mechanism: " *
                "constraint row $i reduces to " *
                "0 = $(wrhs[i]) * log(Keq)")
        end
        push!(pivot_entries, (i, best_col))
        push!(pivot_col_set, best_col)
        pv = wA[i, best_col]
        wA[i, :] ./= pv
        wrhs[i] /= pv
        for r in 1:nc
            if r != i && wA[r, best_col] != 0
                f = wA[r, best_col]
                wA[r, :] .-= f .* wA[i, :]
                wrhs[r] -= f * wrhs[i]
            end
        end
    end

    dep_exprs = Dict{Symbol, Union{Symbol, Expr}}()
    for (prow, pcol) in pivot_entries
        factors = [
            (columns[c], -wA[prow, c])
            for c in 1:n_vars
            if c != pcol && wA[prow, c] != 0
        ]
        dep_exprs[columns[pcol]] = build_power_expr(wrhs[prow], factors)
    end
    dep_set = Set(keys(dep_exprs))
    return dep_exprs, Tuple(p for p in columns if p ∉ dep_set)
end

"""
Gaussian-elimination kernel underlying `_dependent_param_exprs`: assemble the
constraint system and solve it for the dependent/independent partition. Takes
the rename map as a parameter so callers can supply either the user-defined
kinetic-group rename (Pass 1 only) or the full rename map that also absorbs
single-symbol Wegscheider RE ties (Pass 1 + Pass 2).

Pass 2 of `_build_wegscheider_rename_map` calls this kernel with the
Pass-1-only rename to discover which single-symbol ties to absorb; the
display path in `rate_equation_string` likewise calls it with the
Pass-1-only rename to keep absorbed ties visible under the
`# Wegscheider constraints:` section.

`step_params` and `all_params` default to the mechanism's own `:None`-state
symbols. The allosteric per-state derivation passes state-tagged versions so
`name(p, mech)` renders `K_A_…`/`K_I_…`/bare-`:EqualAI` symbols and `all_params`
carries the matching tagged column set (they must agree symbol-for-symbol).
"""
function _dependent_param_exprs_kernel(
    mech::Mechanism,
    rename::AbstractDict{Symbol, Symbol};
    step_params = _step_parameters(mech),
    all_params = _raw_param_symbols(mech),
)
    A, rhs, columns, priority =
        _assemble_constraints(mech, rename; step_params, all_params)
    return _solve_dependent_set(A, rhs, columns, priority)
end

# Type-dispatching wrapper preserves the existing call sites in
# _dependent_param_exprs and _build_kinetic_rename_map / _build_wegscheider_rename_map.
_dependent_param_exprs_kernel(M::Type{<:EnzymeMechanism},
                              rename::AbstractDict{Symbol, Symbol}) =
    _dependent_param_exprs_kernel(Mechanism(M()), rename)

# ─── Preamble Building Helpers ───────────────────────────────────

"""Build destructuring Expr: (; a, b, c) = source"""
function _destructuring_expr(syms, source::Symbol)
    Expr(:(=), Expr(:tuple, Expr(:parameters, syms...)), source)
end

"""
Collect raw parameter symbols (one K or k_f/k_r per kinetic group) plus
`E_total`, in step order.
"""
function _sorted_raw_param_symbols(M::Type{<:EnzymeMechanism})
    Tuple((_raw_param_symbols(M())..., :E_total))
end

"""Full mode: destructure all params + concs, then raw expr."""
function _build_rate_body(M, ::Type{FullMode})
    expr, all_params, conc_syms = _raw_rate_expr_and_symbols(M)
    Expr(:block,
        _destructuring_expr(all_params, :params),
        _destructuring_expr(conc_syms, :concs),
        expr)
end

"""Reduced mode: destructure indep params + concs, define dep params, then raw expr."""
function _build_rate_body(M, ::Type{ReducedMode})
    expr, _, conc_syms = _raw_rate_expr_and_symbols(M)
    dep_exprs, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(conc_syms, :concs),
        assignments...,
        expr)
end
