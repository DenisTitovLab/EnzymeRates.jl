# ABOUTME: Cycle-space helpers for the :OnlyA guard and the split move's counter:
# ABOUTME: the exact null space, the cycle basis and the parameter count of a regrouping.

# ─── Thermodynamic Constraint Infrastructure ─────────────────────

"""
Reduced row echelon form over `Rational{BigInt}`. Returns the pivot and
free column indices (pivot_cols in row-pivot order, so pivot_cols[i] is the
pivot at reduced-matrix row i) plus the reduced matrix R. Used by
`_rational_nullspace` (nullspace basis) and `_partition_independent_count` (rank).
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

_thermodynamic_constraints(mech::Mechanism) =
    _thermodynamic_constraints(mech.reaction, steps(mech))

function _thermodynamic_constraints(rxn::EnzymeReaction, groups::Vector{Vector{Step}})
    flat = _flat_steps(groups)
    ras = reactants(rxn)

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
        c = nu[j0] // nu_net[j0]
        nu == c .* nu_net && isinteger(c) ||
            error("Cycle $i produces metabolite change not proportional to net reaction")
        Int(c)
    end
    return C, Int[classify(i) for i in axes(C, 1)]
end

"""How a step's constants enter the constraint columns: `:ss` (a forward and a
reverse rate), `:binding_K` (a dissociation constant, whose column carries a sign
flip) or `:iso_K` (an equilibrium constant)."""
_count_kind(s::Step) = is_equilibrium(s) ? (is_binding(s) ? :binding_K : :iso_K) : :ss

"""
Add flat step `j`'s cycle incidence `C[:, j]` into the constraint matrix `A` at the
columns `cols` of its constants, by the step's constant kind (`_count_kind`). A
binding K enters negated, because it is a Kd in the polynomial while the cycle
product uses 1/Kd; an iso K enters as is; an SS step contributes `+kf` at `cols[1]`
and `-kr` at `cols[2]`.
"""
function _add_step_column!(A, C, j, kind::Symbol, cols)
    for i in axes(C, 1)
        c = C[i, j]
        c == 0 && continue
        if kind === :ss
            A[i, cols[1]] += c
            A[i, cols[2]] -= c
        else
            A[i, cols[1]] += kind === :binding_K ? -c : c
        end
    end
end

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
constructed child: the joint solve's independent set is its columns of the law less
its pivots, and a kinetic group's constants form one column here and one class there.
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
        for (j, g) in enumerate(group_of_step)
            kind = kind_of_step[j]
            cols = [column[(g, k)] for k in 1:(kind === :ss ? 2 : 1)]
            _add_step_column!(A, C, j, kind, cols)
        end
        length(column) - length(_rref_partition(A)[1])
    end
    counter
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
graph. The check graph is read for its cycles only, so the kept groups enter as
stored, not as a `Mechanism`: dropping a group that mixes a fused binding with
abortive bindings can leave abortive complexes in a rapid-equilibrium segment with
no bottom form, which the `Mechanism` constructor rejects (`_bottomless_re_segment`).
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

Names the offending bindings through the plain `Mechanism` of `cat_steps`; it must
not construct an `AllostericMechanism`, whose constructor calls it.
"""
function _onlya_haldane_violation(rxn::EnzymeReaction,
                                  cat_steps::Vector{Vector{Step}},
                                  cat_allo_states::Vector{Symbol})
    onlya(g) = cat_allo_states[g] === :OnlyA
    keep = [g for g in eachindex(cat_steps)
            if !(onlya(g) && any(_is_chemistry, cat_steps[g]))]
    any(onlya, keep) || return nothing
    kept = cat_steps[keep]
    C, _ = _thermodynamic_constraints(rxn, kept)
    flat = _flat_steps(kept)
    groups = [k for (k, g) in enumerate(keep) if onlya(g)]
    # Column k holds `:OnlyA` group k's `ε` exponents: its steps' cycle entries summed.
    M = Rational{BigInt}[sum((C[i, j] for (j, (_, g)) in enumerate(flat) if g == k);
                             init = 0) for i in axes(C, 1), k in groups]
    cm = Mechanism(rxn, cat_steps)
    labels = [string(name(first(_step_constants(first(kept[k]), :None)), cm))
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
