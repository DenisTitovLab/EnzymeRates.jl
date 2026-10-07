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

_raw_param_symbols(@nospecialize(m::EnzymeMechanism)) = _raw_param_symbols(Mechanism(m))

"""
For each step in `m` (in flat-iteration order), yield the Parameter
instances that govern that step: `[Kd|Kiso]` for an RE step, `[Kon|Kfor,
Koff|Krev]` for an SS step. Each Parameter is anchored on the original
step (not the rep), so `name(p, m)` renders to the rep's structural
Symbol via the value-context chokepoint, collapsing kinetic-group members.
"""
_step_parameters(m::Mechanism) =
    Vector{Parameter}[_step_constants(s, :None) for (s, _) in _flat_steps(m)]

# ─── Structural primacy: free-enzyme set + step priority ─────────

"""
Set of enzyme-form names that are NOT the RHS of any RE step that
consumes a metabolite, `F + met… ⇌ F_bound`. Walks `Mechanism.steps`
directly: such a step leaves its consumed metabolites on `to_species`, so
`to_species`'s name is excluded from the free set. Iso steps don't
determine binding state, and SS steps take no part in this rule (they do in
the next one).

A form that carries bound metabolites but has no consuming step into it in this
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
function _compute_free_enz_set(m::Union{Mechanism, AllostericMechanism})
    flat = [s for g in steps(m) for s in g]
    # A step that consumes a metabolite leaves it on to_species. The from-side is
    # the "free + met" reactant; the to-side is the bound form.
    into(re) = Set{Symbol}(name(to_species(s)) for s in flat
                           if !isempty(consumed(s)) && (!re || is_equilibrium(s)))
    bound_in, re_bound_in = into(false), into(true)
    Set{Symbol}(name(sp) for s in flat for sp in (from_species(s), to_species(s))
                if name(sp) ∉ re_bound_in && (isempty(bound(sp)) || name(sp) in bound_in))
end

"""
The free-enzyme form names of `m` (`_compute_free_enz_set`), cached per mechanism.
Callers share the returned set and must not mutate it.
"""
function _free_enz_set(m::Union{Mechanism, AllostericMechanism})
    c = m.naming
    fes = c.free_enz
    fes === nothing || return fes
    c.free_enz = _compute_free_enz_set(m)
end

"""
Structural primacy base score for a step (lower = more primary / less
eliminable). A metabolite step takes up or gives off a free metabolite (a pure
binding, a fused step or a Theorell–Chance step). Free-enzyme RE metabolite
step (-1) < free-enzyme SS metabolite step (0) < non-free metabolite step
(10) < internal isomerization (20). Shared by the
kinetic-group name representative (argmin) and the Haldane elimination pivot
(argmax, which adds a +0/+1 forward/reverse offset per rate constant).
"""
function _step_priority(s::Step, free_enz_set::Set{Symbol})
    has_met = !is_iso(s)
    is_free = (name(from_species(s)) in free_enz_set) ||
              (name(to_species(s))   in free_enz_set)
    is_equilibrium(s) && has_met && is_free && return -1
    return !has_met ? 20 : is_free ? 0 : 10
end

"""
Kinetic-group naming representative: the structurally-primary step
(`argmin _step_priority`), with a deterministic lexical tiebreak between two
distinct steps (`_step_canonical_key`: species pair + consumed and released
metabolites + RE/SS flag).
"""
_group_rep(group::Vector{Step}, free_enz_set::Set{Symbol}) =
    argmin(s -> (_step_priority(s, free_enz_set), _step_canonical_key(s)), group)

# ─── Thermodynamic Constraint Infrastructure ─────────────────────

"""
Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
pivot at reduced-matrix row i) plus the reduced matrix R. Used by
`_rational_nullspace` (nullspace basis), `_partition_independent_count` (rank) and
`_solve_dependent_set` (dependent parameters).
"""
function _rref_partition(A::AbstractMatrix)
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

"""
Basis of `{x : M·x = 0}` as the COLUMNS of the returned matrix, exact over the
rationals: each free column of `M`'s reduced row echelon form yields one basis
vector, 1 at its free coordinate.
"""
function _rational_nullspace(M::AbstractMatrix)
    pivots, free, R = _rref_partition(M)
    N = zeros(Rational{BigInt}, size(M, 2), length(free))
    for (k, f) in enumerate(free)
        N[f, k] = 1
        for (r, c) in enumerate(pivots); N[c, k] = -R[r, f]; end
    end
    N
end

"""
`_rational_nullspace(A)` with each basis column scaled to coprime integers whose
first nonzero entry is positive.
"""
function _integer_nullspace(A::Matrix{Int})
    N = _rational_nullspace(A)
    result = zeros(Int, size(N))
    for k in axes(N, 2)
        col = Int.(N[:, k] .* lcm(denominator.(N[:, k])))
        col .÷= gcd(col)
        result[:, k] .= col[findfirst(!iszero, col)] < 0 ? -col : col
    end
    result
end

function _thermodynamic_constraints(mech::Mechanism)
    flat = _flat_steps(mech)
    ras = reactants(mech.reaction)

    # Enzyme incidence matrix (rows = enzyme forms in first-seen step-walk order)
    form = Dict{Symbol, Int}()
    for (s, _) in flat, sp in (from_species(s), to_species(s))
        get!(form, name(sp), length(form) + 1)
    end
    B = zeros(Int, length(form), length(flat))
    for (j, (s, _)) in enumerate(flat)
        B[form[name(from_species(s))], j] -= 1
        B[form[name(to_species(s))], j] += 1
    end

    # Stoichiometry matrix (rows = metabolites, cols = steps), read from each
    # step's consumed (-1) and released (+1) lists.
    #
    # Iso steps carry no free-pool metabolite — their bound content is
    # encoded in the enzyme-form identity — so both lists are empty and they
    # contribute zero. Do NOT add a from_bound/to_bound diff for iso steps:
    # that double-counts metabolites already accounted for by the
    # binding/release steps and inflates the cycle's net change (e.g.
    # 1/Keq -> 1/Keq^2).
    met_idx = Dict(name(metabolite(ra)) => i for (i, ra) in enumerate(ras))
    stoich_mat = zeros(Int, length(ras), length(flat))
    for (j, (s, _)) in enumerate(flat), (mets, sgn) in ((consumed(s), -1), (released(s), 1))
        for m in mets
            haskey(met_idx, name(m)) && (stoich_mat[met_idx[name(m)], j] += sgn)
        end
    end
    nu_net = zeros(Int, length(ras))
    for ra in ras
        nu_net[met_idx[name(metabolite(ra))]] += metabolite(ra) isa Substrate ? -1 : 1
    end

    # Classify each null-space cycle as Haldane (proportional to the
    # net reaction → contributes log(Keq)) or Wegscheider (closed
    # cycle, zero net change). Errors on cycles that touch metabolites
    # but aren't proportional to the net reaction.
    C = _integer_nullspace(B)'
    j0 = findfirst(!iszero, nu_net)
    function classify(i)
        nu = stoich_mat * C[i, :]
        c = j0 === nothing ? 0 // 1 : nu[j0] // nu_net[j0]
        nu == c .* nu_net && isinteger(c) ||
            error("Cycle $i produces metabolite change not proportional to net reaction")
        Int(c)
    end
    return C, Int[classify(i) for i in axes(C, 1)]
end

"""
    _dependent_param_exprs(mech::Mechanism) → (dep_exprs, indep_params)

Select dependent parameters and build substitution expressions for the
Haldane / Wegscheider thermodynamic constraints. Steps in the same
kinetic group share parameters: their cycle-incidence columns are merged
into the representative step's column before Gaussian elimination, so
`dep_exprs` and `indep_params` are keyed only on representatives.

Calls `_build_wegscheider_rename_map(mech)` to obtain the rename map for
absorbed single-symbol Wegscheider RE ties and forwards to
`_dependent_param_exprs_kernel`. The `Type{<:AbstractEnzymeMechanism}` method lifts
with `_concrete` and delegates here.
"""
function _dependent_param_exprs(mech::Mechanism)
    rename = _build_wegscheider_rename_map(mech)
    dep_exprs, indep = _dependent_param_exprs_kernel(mech, rename)
    # Filter Pass-2-absorbed symbols out of indep. Pass 2 of
    # `_build_wegscheider_rename_map` adds entries like `K_EP_to_E_P => K_ES_to_E_S`
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

"""Number of independent (fitted) rate constants of a concrete mechanism, computed
from the thermodynamic constraint solve without compiling the mechanism."""
_independent_param_count(m::Union{Mechanism, AllostericMechanism}) =
    length(_dependent_param_exprs(m)[2])

"""How a step's constants enter the constraint columns: `:ss` (a forward and a
reverse rate), `:binding_K` (a dissociation constant, whose column carries a sign
flip) or `:iso_K` (an equilibrium constant)."""
_count_kind(s::Step) = is_equilibrium(s) ? (is_binding(s) ? :binding_K : :iso_K) : :ss

"""
    _partition_independent_count(parent::Mechanism) -> counter

Return `counter(group_of_step, kind_of_step = kinds of the parent's steps)`, the
independent-parameter count of the mechanism obtained by regrouping `parent`'s
flat steps (in `_flat_steps` order) into the groups labelled by `group_of_step`,
each step counted under the constant kind `kind_of_step` gives it (`_count_kind`).
Regrouping moves no edges and a flag changes none, so the cycle basis of the step
graph (`_thermodynamic_constraints`) is the same for every call and is computed
once here; each call only merges step columns by group and takes the rank. The
split move passes a step's own kind, or `:binding_K`/`:iso_K` for a steady-state
step it reverts to rapid equilibrium. Equals `_independent_param_count` of the
constructed child: the kernel's independent set is the columns minus the pivots,
and folding a single-symbol Wegscheider tie onto its target removes one column and
one rank together, so the count is invariant to the rename. Column sign
conventions match `_assemble_constraints`: a binding K enters with a sign flip, an
iso K without, an SS step contributes `+kf` and `-kr`.
"""
function _partition_independent_count(parent::Mechanism)
    C, _ = _thermodynamic_constraints(parent)
    kinds = [_count_kind(s) for (s, _) in _flat_steps(parent)]
    function counter(group_of_step::AbstractVector{Int},
                     kind_of_step::AbstractVector{Symbol} = kinds)
        length(group_of_step) == length(kinds) ||
            error("group_of_step must label every flat step of the parent")
        length(kind_of_step) == length(kinds) ||
            error("kind_of_step must label every flat step of the parent")
        column = Dict{Tuple{Int, Int}, Int}()
        for (j, g) in enumerate(group_of_step)
            get!(column, (g, 1), length(column) + 1)
            kind_of_step[j] === :ss && get!(column, (g, 2), length(column) + 1)
        end
        A = zeros(Int, size(C, 1), length(column))
        for (j, g) in enumerate(group_of_step), i in axes(C, 1)
            c = C[i, j]
            c == 0 && continue
            if kind_of_step[j] === :ss
                A[i, column[(g, 1)]] += c
                A[i, column[(g, 2)]] -= c
            else
                A[i, column[(g, 1)]] += kind_of_step[j] === :binding_K ? -c : c
            end
        end
        length(column) - length(_rref_partition(A)[1])
    end
    counter
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

    step_name(p::Parameter) = (s = name(p, mech); get(rename, s, s))

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
    # No value is a never-pivot sentinel.
    priority = fill((is_i_state, 0), n_vars)
    for j in 1:nsteps
        step = step_params[j][1].step
        base = _step_priority(step, free_enz_set)
        if is_equilibrium(flat[j][1])
            s = step_name(step_params[j][1])
            haskey(sym_col, s) && (priority[sym_col[s]] = (is_i_state, base))
        else
            # A steady-state step's reverse constant is eliminated before its forward one.
            # A fused binding of a product is a release stored as the binding it reverses,
            # so its forward constant runs against the reaction and goes first instead,
            # keeping the catalytic constant fitted.
            against = _is_chemistry(step) && bound_metabolite(step) isa Product
            for (offset, p) in enumerate(step_params[j])
                s = step_name(p)
                rank = against ? 2 - offset : offset - 1
                haskey(sym_col, s) && (priority[sym_col[s]] = (is_i_state, base + rank))
            end
        end
    end
    return A, rhs, columns, priority
end

"""
True when some `y` makes every row of `N·y` strictly positive — i.e. the open
cone `{y : N·y > 0}` is nonempty. Exact Fourier-Motzkin elimination: to drop
variable `v`, every (positive, negative) row pair is combined with positive
coefficients, which preserves strictness; a row that reduces to all-zero
encodes `0 > 0` and refutes the system. Returns `nothing` if elimination blows
up, so the caller can fall back to the sound per-row test rather than error on
the enumeration hot path.
"""
function _has_strict_positive_combination(N::AbstractMatrix{Rational{BigInt}})
    d = size(N, 2)
    d == 0 && return false
    rows = [collect(N[i, :]) for i in axes(N, 1)]
    for v in d:-1:1
        pos = [r for r in rows if r[v] > 0]
        neg = [r for r in rows if r[v] < 0]
        nxt = [r for r in rows if r[v] == 0]
        length(nxt) + length(pos) * length(neg) > 4000 && return nothing
        for p in pos, q in neg
            push!(nxt, p .* (-q[v]) .+ q .* p[v])
        end
        any(r -> all(iszero, r), nxt) && return false
        rows = nxt
    end
    true
end

_onlya_violation_message(offenders) =
    "an :OnlyA binding ($(join(offenders, ", "))) leaves a " *
    "thermodynamic (Haldane/Wegscheider) cycle unsatisfiable: the " *
    "inactive conformation cannot close that cycle at finite nonzero " *
    "affinity. Tag the cycle's chemical step :OnlyA, or tag an " *
    "opposing binding :OnlyA so the affinities diverge together."

"""
    _onlya_haldane_violation(rxn, cat_steps, cat_allo_states)
        → Union{Nothing, String}

Return `nothing` when every catalytic thermodynamic cycle's Haldane relation
stays satisfiable under the `K_I = K_A/ε`, `ε → 0⁺` limit an `:OnlyA` binding
asserts; otherwise return a message naming the offending `:OnlyA` bindings.

A cycle's Haldane carries `∏ε_p/∏ε_s` on the inactive side. The `ε` are
independent, so that monomial can be held finite only when its exponents
carry both signs. All-same-sign drives it to `0` or `∞`; the only thing that
absorbs the imbalance is a free inactive catalytic ratio `k_I_f/k_I_r`. An
`:OnlyA` catalytic tag supplies exactly that — its rate constants vanish from
the rate equation, so their ratio is free to satisfy the cycle's Haldane at any
affinity — which is why the inactive Haldane is present but never binding here.

The check graph drops `:OnlyA` chemistry groups (those holding a chemistry step,
`_is_chemistry`: an isomerization, a fused binding or a Theorell–Chance step), so
a cycle running through one never appears and never reports a violation: that
free `k_I` ratio is the escape. The `:OnlyA` bindings are the `:OnlyA` groups of
plain bindings alone; a fused binding tagged `:OnlyA` is chemistry and leaves the
graph.
Bindings completing no cycle (competitive inhibitors, dead ends, regulator
sites) never enter a row and take no part. Both catalytic (Haldane) and
binding-only (Wegscheider, `rhs = 0`) cycle rows are inspected, so a one-sided
`:OnlyA` binding on a pure random-order binding square is caught.

The verdict comes from two stages. The per-row sign test runs first: it flags a
cycle whose single row's `:OnlyA` exponents are all one sign, which forces a sum
of same-signed positive terms to vanish and so is a sound rejection. It is cheap
and keeps the common rejections on the fast path, but it reads one row at a time
and a multi-cycle coupled inconsistency passes it. The complete condition
follows: writing `ε_p = δ^{w_p}`, every cycle's monomial is finite nonzero for
some `δ → 0⁺` exactly when the coupled `ε`-exponent system admits a
strictly-positive nullspace vector, `M·w = 0` with `w > 0` (Stiemke
feasibility). The two stages agree for every random-order mechanism up to bi-bi;
they part from ter-substrate up.

An `:OnlyA` group's `ε` exponent in a cycle is the sum of its steps' entries in
that row of the cycle basis (`_thermodynamic_constraints`), for RE and SS
bindings alike: an RE binding's `1/Kd` and an SS binding's `kon/koff` both enter
the cycle's product raised to the step's entry, so a cycle that mixes the two step
kinds reads its signs on one scale and a balanced pair is accepted.

Builds a plain `Mechanism`; it must not call `_state_allo_mechanism`, which
would construct an `AllostericMechanism` and recurse.
"""
function _onlya_haldane_violation(rxn::EnzymeReaction,
                                  cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    onlya(g) = cat_allo_states[g] === :OnlyA
    keep = [g for g in eachindex(cat_steps)
            if !(onlya(g) && any(_is_chemistry, cat_steps[g]))]
    any(onlya, keep) || return nothing
    onlyA_steps = Set{Step}(s for g in keep if onlya(g) for s in cat_steps[g])
    cm = Mechanism(rxn, [copy(cat_steps[g]) for g in keep])
    C, _ = _thermodynamic_constraints(cm)
    flat = _flat_steps(cm)
    groups = unique(g for (s, g) in flat if s in onlyA_steps)
    isempty(groups) && return nothing
    # Column k holds `:OnlyA` group k's `ε` exponents: its steps' cycle entries summed.
    M = Rational{BigInt}[sum((C[i, j] for (j, (_, g)) in enumerate(flat) if g == k);
                             init = 0) for i in axes(C, 1), k in groups]
    labels = [string(name(first(_step_constants(first(steps(cm)[k]), :None)), cm))
              for k in groups]
    # Per-row sign test first: sound (an all-one-sign row forces a sum of
    # same-signed positive terms to vanish) and cheap, so it keeps the common
    # rejections on the fast path.
    for row in eachrow(M)
        nz = findall(!iszero, row)
        !isempty(nz) && allequal(sign.(row[nz])) &&
            return _onlya_violation_message(sort!(labels[nz]))
    end

    # The complete condition: the `ε` exponents must admit a strictly-positive
    # solution of `M·w = 0` (Stiemke feasibility). The per-row test above only
    # inspects one row at a time, so from ter-substrate up a multi-cycle coupled
    # inconsistency passes it. `nothing` from the cone test means elimination
    # blew up; fall back to the per-row verdict (sound, just incomplete) rather
    # than reject a mechanism we could not decide.
    _has_strict_positive_combination(_rational_nullspace(M)) === false || return nothing
    return _onlya_violation_message(sort(labels))
end

"""
Solve an assembled constraint system `(A, rhs, columns, priority)` for a
dependent/independent parameter partition. The columns are sorted by priority,
highest first (every `is_i_state` column ahead of the rest, ties in ascending column
order), and `[A rhs]` is brought to reduced row echelon form. Each pivot column
becomes dependent, expressed via the non-pivot columns of its row; every other
column is independent (fitted). A pivot in the `rhs` column means the rows combine
to `0 = log Keq`, a thermodynamic contradiction, which errors. Returns
`(dep_exprs, indep)`. An empty system (no rows) yields no dependents and all
columns independent.

The pivots are the greedy basis over the sorted columns: a column becomes dependent
exactly when it is not a linear combination of the higher-priority columns, so the
constraints eliminate the highest-priority columns they can. Eliminating row by row
and pivoting each row on its highest-priority remaining column picks the same set,
and the reduced row echelon form is unique, so the dependent expressions match too.
"""
function _solve_dependent_set(
    A::AbstractMatrix{Rational{BigInt}},
    rhs::AbstractVector{Rational{BigInt}},
    columns::AbstractVector{Symbol},
    priority::AbstractVector{Tuple{Bool, Int}},
)
    n_vars = length(columns)
    order = sortperm(priority; rev = true)
    pivots, _, R = _rref_partition([A[:, order] rhs])
    n_vars + 1 in pivots && error(
        "Thermodynamically contradictory mechanism: a combination of its " *
        "constraint rows reduces to 0 = log(Keq)")
    dep_exprs = Dict{Symbol, Union{Symbol, Expr}}()
    for (r, p) in enumerate(pivots)
        factors = [(columns[order[c]], -R[r, c])
                   for c in 1:n_vars if c != p && R[r, c] != 0]
        dep_exprs[columns[order[p]]] = build_power_expr(R[r, n_vars + 1], factors)
    end
    return dep_exprs, Tuple(p for p in columns if !haskey(dep_exprs, p))
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

"""
Type-dispatching wrapper preserves the existing call sites in
_dependent_param_exprs and _build_kinetic_rename_map / _build_wegscheider_rename_map.
"""
_dependent_param_exprs_kernel(@nospecialize(M::Type{<:EnzymeMechanism}),
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
function _sorted_raw_param_symbols(@nospecialize(M::Type{<:EnzymeMechanism}))
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
