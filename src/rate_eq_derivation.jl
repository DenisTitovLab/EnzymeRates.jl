# ABOUTME: King-Altman/Cha rate equation derivation via @generated functions.
# ABOUTME: Parameters API, kcat computation, and MWC allosteric rate assembly.

# ─── Parameters API ─────────────────────────────────────────

const _AnyMechanism = AbstractEnzymeMechanism

"""
Suffix appended to single-symbol equality lines whose LHS got folded
into the kinetic-group rename map. Both display sites (User defined
kinetic-group merges and absorbed single-symbol Wegscheider ties) emit
this exact string so the rate-equation dedup key in
`identify_rate_equation.jl` can strip these provenance lines.
"""
const ANNOTATION_SUBSTITUTED = "  (substituted into v)"

"""
    parameters(m::EnzymeMechanism, [mode])
    parameters(m::AllostericEnzymeMechanism, [mode])

Return the parameter names required for the given mode as a tuple of Symbols.

# Modes
- `Reduced` (default): independent k's + Keq + E_total. The set of
  symbols the user supplies to evaluate the Haldane-reduced rate
  equation. Returned for both `EnzymeMechanism` and
  `AllostericEnzymeMechanism`.
- `Full`: all raw rate-constant symbols + `E_total`. For
  `EnzymeMechanism` this is "all 2N k's + `E_total`." For
  `AllostericEnzymeMechanism` it composes the catalytic raw A-state
  symbols + every I-state mirror (catalytic + regulatory + synthesized
  dep) + reg-site A-state K's (skipping `:OnlyI` ligands) + `:L` +
  `:E_total`. The allosteric Full mode enumerates the complete raw-symbol
  set; no `rate_equation` method is defined for
  `(::AllostericEnzymeMechanism, ::FullMode)`, so this mode is for
  symbol enumeration, not runtime evaluation.
"""
function parameters end

parameters(m::_AnyMechanism) = parameters(m, Reduced)

# ── EnzymeMechanism ───────────────────────────────────────────
@generated function parameters(
    ::EnzymeMechanism{Sig}, ::FullMode,
) where {Sig}
    mech = Mechanism(EnzymeMechanism{Sig}())
    params = _enumerate_parameters_full(mech)
    names = Tuple(name(p, mech) for p in params)
    Tuple((names..., :E_total))
end

@generated function parameters(::M, ::ReducedMode) where {M <: EnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

# ── AllostericEnzymeMechanism ────────────────────────────────
@generated function parameters(
    ::M, ::FullMode,
) where {M <: AllostericEnzymeMechanism}
    aem = M()
    am = AllostericMechanism(aem)
    # The full symbol set over-emits an I-state mirror for every non-`:OnlyA`
    # catalytic group (`_all_i_state_parameters`). A Case-B synthesized
    # dependent — an `:EqualAI` dep whose Haldane RHS references a
    # `:NonequalAI` symbol, e.g. PK's `k_I_EATPPyruvate_to_EADPPEP` — is the
    # I-form of an `:EqualAI` catalytic reverse rate, so it always coincides
    # with that group's over-emitted `Kfor`/`Krev`(:I) mirror already in
    # `names`. No separate synthesized-name splice is needed.
    params = _enumerate_parameters_full_allosteric(am)
    names = Symbol[name(p, am) for p in params]
    Tuple((names..., :E_total))
end

@generated function parameters(
    ::M, ::ReducedMode,
) where {M <: AllostericEnzymeMechanism}
    _, indep = _dependent_param_exprs(M)
    (indep..., :Keq, :E_total)
end

"""Independent rate constant names for fitting (excludes Keq, E_total)."""
@generated function fitted_params(::M) where {M <: _AnyMechanism}
    _, indep = _dependent_param_exprs(M)
    indep
end

# Convenience methods for the concrete Mechanism / AllostericMechanism
# working representation: lift to the compiled singleton via
# `compile_mechanism`, then dispatch to the singleton method. Callers that
# hold a `Mechanism` (e.g. the `identify_rate_equation` pipeline) use these
# without an explicit lift. `compile_mechanism` does Sig conversion (it
# allocates) — it is NOT the rate_equation hot path. Pass an
# `EnzymeMechanism` / `AllostericEnzymeMechanism` for the 0-alloc/<100ns
# `rate_equation` guarantee.
parameters(m::Union{Mechanism, AllostericMechanism},
           mode::AbstractRateEquationMode) =
    parameters(compile_mechanism(m), mode)
parameters(m::Union{Mechanism, AllostericMechanism}) = parameters(m, Reduced)
fitted_params(m::Union{Mechanism, AllostericMechanism}) =
    fitted_params(compile_mechanism(m))

"""
Build a renaming map for single-symbol Wegscheider RE ties between two
binding K's. Calls `_dependent_param_exprs_kernel` to discover
binding-K-to-binding-K Wegscheider closures of the form `K_a = K_b`
(RHS is a bare Symbol). Both sides must be binding K's (RE step with
metabolite on LHS) — absorbing a binding-K-to-iso-K tie would produce
inconsistent sign-flips when the kernel runs with the full rename,
since the binding-K column is sign-flipped (Kd convention) but the
iso-K column is not.

The rename means the polynomial in `v` uses the representative symbol
directly, so Source-C duplicates (split kinetic groups that
Wegscheider ties back together) collapse at hash time.
"""
function _build_wegscheider_rename_map(mech::Mechanism)
    rename = Dict{Symbol, Symbol}()
    step_params = _step_parameters(mech)
    # binding-K set: value-context rep name of each RE binding step. Walk
    # Mechanism.steps directly — an RE step carrying a bound metabolite is
    # a binding step; step_params is indexed in the same flat order.
    binding_set = Set{Symbol}()
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_set, name(step_params[idx][1], mech))
    end
    # Pass 2: single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _dependent_param_exprs_kernel(mech, rename)
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

_build_wegscheider_rename_map(M::Type{<:EnzymeMechanism}) =
    _build_wegscheider_rename_map(Mechanism(M()))
_build_wegscheider_rename_map(m::EnzymeMechanism) =
    _build_wegscheider_rename_map(typeof(m))

# ─── RE Group Helpers ───────────────────────────────────────

"""
Per-step side breakdown for the rate equation derivation. Returns
`(from_species_sym, to_species_sym, m_lhs_syms, m_rhs_syms)` for a
single `Step`. This is the canonical metabolite-on-which-side
projection: it reads Step fields directly, placing the bound
metabolite on the binding (m_lhs) or release (m_rhs) side from the
canonical metabolite-on-`to_species` placement plus the SS-dissociation
rule. The five-branch logic is load-bearing: SS catalytic-release steps
where the bound metabolite is a Product that doesn't appear in either
bound list put the metabolite on m_rhs, not m_lhs.
"""
function _step_sides(s::Step)
    e_lhs = name(from_species(s))
    e_rhs = name(to_species(s))
    is_iso(s) && return (e_lhs, e_rhs, Symbol[], Symbol[])
    bm = bound_metabolite(s)
    bm_name = name(bm)
    from_bound_names = Set{Symbol}(name(m) for m in bound(from_species(s)))
    to_bound_names = Set{Symbol}(name(m) for m in bound(to_species(s)))
    if bm_name in to_bound_names
        # Canonical binding: bound met on to_species → emit on m_lhs
        return (e_lhs, e_rhs, Symbol[bm_name], Symbol[])
    elseif bm_name in from_bound_names
        # Reverse-canonical (defensive: ctor swap should have prevented this)
        return (e_lhs, e_rhs, Symbol[], Symbol[bm_name])
    elseif !is_equilibrium(s) && bm isa Product &&
           !(isempty(from_bound_names) && isempty(to_bound_names))
        # SS dissociation rule: bound_metabolite is a Product released
        # in an SS catalytic step; not in either bound list (because the
        # Species canonicalization moved it). Emit on m_rhs.
        return (e_lhs, e_rhs, Symbol[], Symbol[bm_name])
    elseif length(from_bound_names) > length(to_bound_names)
        # Bound-list-size fallback: the side with fewer bound metabolites
        # is the release side, so the metabolite goes on m_rhs.
        return (e_lhs, e_rhs, Symbol[], Symbol[bm_name])
    else
        return (e_lhs, e_rhs, Symbol[bm_name], Symbol[])
    end
end

"""
Compute RE-connected groups via union-find over enzyme Species. Walks
`m.steps` directly; nodes are Species values (compared by `==`). Returns
`(enz_species, groups, form_to_group)` with `groups[g] :: Vector{Int}`
indexing into `enz_species`.
"""
function _compute_re_groups(mech::Mechanism)
    enz_species = _enumerate_species(mech)
    N = length(enz_species)
    parent = collect(1:N)
    function find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end
        x
    end
    for group in steps(mech)
        for s in group
            is_equilibrium(s) || continue
            i_from = findfirst(==(from_species(s)), enz_species)
            i_to   = findfirst(==(to_species(s)),   enz_species)
            ra, rb = find(i_from), find(i_to)
            ra != rb && (parent[ra] = rb)
        end
    end
    root_to_group = Dict{Int, Int}()
    groups = Vector{Vector{Int}}()
    form_to_group = zeros(Int, N)
    for i in 1:N
        r = find(i)
        g = get!(root_to_group, r) do; push!(groups, Int[]); length(groups) end
        push!(groups[g], i); form_to_group[i] = g
    end
    enz_species, groups, form_to_group
end

"""Distinct enzyme Species in `m.steps`, in step-walk order."""
function _enumerate_species(m::Mechanism)
    seen = Species[]
    for group in steps(m)
        for s in group
            from_species(s) in seen || push!(seen, from_species(s))
            to_species(s)   in seen || push!(seen, to_species(s))
        end
    end
    seen
end

"""
Concentration symbols (substrates ∪ products ∪ regulators) that may appear in
the rate-equation polynomials — the same set `_poly_to_expr` treats as
concentrations and the only symbols `_reduce_conc_lowest_terms` is allowed to
shift. Parameter symbols are never in this set, so the concentration-GCD can
never drop a fitted parameter.
"""
function _concentration_symbols(mech::Mechanism)
    rxn = reaction(mech)
    cs = Set{Symbol}(name(s) for s in substrates(rxn))
    for p in products(rxn);    push!(cs, name(p)); end
    for rm in regulators(rxn); push!(cs, name(regulator(rm))); end
    cs
end

"""
Free enzyme of a rapid-equilibrium segment: the form with the fewest bound
metabolites, tie-broken toward no covalent residual, then a deterministic
name. Referencing each segment's alphas to this form makes the derivation
independent of the step order and yields the readable `1 + [S]/K + …` form.
"""
function _segment_root(group, enz_species)
    argmin(i -> (length(bound(enz_species[i])),
                 has_residual(enz_species[i]) ? 1 : 0,
                 string(name(enz_species[i]))), group)
end

"""
Compute alpha factors (relative concentrations within RE groups) as POLY
values. Iterates `mech.steps` directly. Binding steps' direction comes
from the canonical Step form (metabolite-on-`to_species`); iso steps
are physical-forward (canonicalized in the Mechanism constructor).
`step_to_K[idx]` is the parameter Symbol for the RE step at flat
position `idx` (rep-renamed via the `name(p, m)` chokepoint).
"""
function _compute_alpha(mech::Mechanism, enz_species,
                        enz_name_to_form, groups, step_to_K)
    N = length(enz_species)
    alpha = Vector{POLY}(fill(poly_one(), N))
    flat = _flat_steps(mech)

    for group in groups
        length(group) == 1 && continue
        root = _segment_root(group, enz_species)
        visited = Set{Int}([root])
        queue = [root]
        while !isempty(queue)
            cur = popfirst!(queue)
            for (idx, (s, _)) in enumerate(flat)
                is_equilibrium(s) || continue
                e_l, e_r, m_l, m_r = _step_sides(s)
                i_f = enz_name_to_form[e_l]
                j_f = enz_name_to_form[e_r]
                Ksym = step_to_K[idx]
                Kp = poly_sym(Ksym)
                Kinv = POLY(_mono(Ksym => -1) => 1)
                is_iso = isempty(m_l) && isempty(m_r)
                if i_f == cur && j_f ∉ visited
                    if is_iso
                        alpha[j_f] = poly_mul(alpha[cur], Kp)
                    else
                        f = poly_mul(poly_sym(m_l[1]), Kinv)
                        isempty(m_r) ||
                            (f = poly_mul(f, POLY(_mono(m_r[1] => -1) => 1)))
                        alpha[j_f] = poly_mul(alpha[cur], f)
                    end
                    push!(visited, j_f); push!(queue, j_f)
                elseif j_f == cur && i_f ∉ visited
                    if is_iso
                        alpha[i_f] = poly_mul(alpha[cur], Kinv)
                    else
                        f = poly_mul(Kp, POLY(_mono(m_l[1] => -1) => 1))
                        alpha[i_f] = poly_mul(alpha[cur], f)
                    end
                    push!(visited, i_f); push!(queue, i_f)
                end
            end
        end
    end
    alpha
end

"""Build rate poly for one SS step direction in Laurent (fractional) form."""
function _ss_contrib(k_poly, mets, i_form, alpha)
    r = isempty(mets) ? k_poly : poly_mul(k_poly, reduce(poly_mul, poly_sym.(mets)))
    poly_mul(r, alpha[i_form])
end

# ─── Raw Rate Equation Derivation (Unified Cha / King-Altman) ───

"""
Build raw numerator and denominator POLYs for the rate equation by
walking the lifted `Mechanism`. Parameter Symbols on the leaves of
`num`/`den` are produced via the `name(p, mech)` chokepoint (which
collapses kinetic-group members to their rep's name). `rename_map` then
applies any single-symbol Wegscheider ties as a post-pass.
"""
function _raw_symbolic_rate_polys(mech::Mechanism, step_params, rename_map,
                                  subs_species, prods_species;
                                  allow_dead::Bool=false)
    enz_species, groups, form_to_group = _compute_re_groups(mech)
    enz_name_to_form = Dict{Symbol, Int}(
        name(es) => i for (i, es) in enumerate(enz_species))
    flat = _flat_steps(mech)
    step_to_K = Dict{Int, Symbol}(
        i => name(step_params[i][1], mech)
        for i in eachindex(flat) if is_equilibrium(flat[i][1]))
    alpha = _compute_alpha(mech, enz_species,
                           enz_name_to_form, groups, step_to_K)
    G = length(groups)

    R = [poly_zero() for _ in 1:G, _ in 1:G]
    for (idx, (s, _)) in enumerate(flat)
        is_equilibrium(s) && continue
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(s)
        i_form = enz_name_to_form[e_lhs]
        j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        kf_poly = poly_sym(name(step_params[idx][1], mech))
        kr_poly = poly_sym(name(step_params[idx][2], mech))
        R[g1, g2] = poly_add(R[g1, g2],
            _ss_contrib(kf_poly, m_lhs, i_form, alpha))
        R[g2, g1] = poly_add(R[g2, g1],
            _ss_contrib(kr_poly, m_rhs, j_form, alpha))
    end

    L = [i == j ? poly_zero() : poly_neg(R[i,j])
         for i in 1:G, j in 1:G]
    for i in 1:G
        L[i, i] = reduce(poly_add, R[i, j] for j in 1:G if j != i; init=poly_zero())
    end
    D = [begin
        idx = [r for r in 1:G if r != root]
        isempty(idx) ? poly_one() : sym_det(L[idx, idx], G - 1)
    end for root in 1:G]

    den = poly_zero()
    for g in 1:G
        sigma = reduce(poly_add, (alpha[i] for i in groups[g]); init=poly_zero())
        csigma = _rename_symbols(sigma, rename_map)
        den = poly_add(den, poly_mul(csigma, D[g]))
    end

    num = _compute_numerator(
        mech, enz_name_to_form, step_params,
        alpha, form_to_group,
        D, subs_species, prods_species; allow_dead=allow_dead)

    num = _rename_symbols(num, rename_map)
    den = _rename_symbols(den, rename_map)
    conc_set = _concentration_symbols(mech)
    num, den = _reduce_conc_lowest_terms(num, den, conc_set)
    num, den
end

function _raw_symbolic_rate_polys(M::Type{<:EnzymeMechanism})
    mech = Mechanism(M())
    step_params = _step_parameters(mech)
    rename_map = _build_wegscheider_rename_map(M)
    # substrates(::EnzymeReaction) returns Vector{Substrate} (concrete metabolite
    # structs); _compute_numerator compares them as Symbols (`name(...) == S`,
    # `S in subs_in(...)`). Wrap explicitly.
    subs_syms = Symbol[name(s) for s in substrates(mech.reaction)]
    prods_syms = Symbol[name(p) for p in products(mech.reaction)]
    _raw_symbolic_rate_polys(mech, step_params, rename_map,
                              subs_syms, prods_syms)
end

"""
Numerator = net flux across one complete steady-state reaction-cut. Each
per-turnover-conserved "event" is a candidate cut whose SS-step fluxes sum to v:
metabolite cuts (bind a substrate / release a product / iso-convert a substrate /
iso-produce a product) and central-species cuts (produce / consume an iso-step
endpoint form). Dead-end binding/release steps touching a substrate-product mixed
complex are excluded (a chemistry step producing such a complex — a ping-pong
covalent intermediate — is kept). A candidate is usable iff all its steps are SS;
NUM = oriented-flux sum over the chosen cut (prefer a metabolite cut — always
complete — over a central-species cut, then fewest steps, then a chemistry cut,
then sorted indices). No usable candidate ⇒ a complete all-RE catalytic cycle ⇒
no finite rate ⇒ raise. A central-species cut is trusted only when its form is
unique up to bound regulators; a regulator-variant sibling is a parallel route the
single-form cut would undercount, so raise rather than return a wrong rate.
"""
function _compute_numerator(
    mech::Mechanism, enz_name_to_form, step_params,
    alpha, form_to_group, D, subs_species, prods_species;
    allow_dead::Bool=false,
)
    is_mixed_substrate_product_complex(f) =
        any(m -> m isa Substrate, bound(f)) && any(m -> m isa Product, bound(f))
    subs_in(f)  = Set(name(m) for m in bound(f) if m isa Substrate)
    prods_in(f) = Set(name(m) for m in bound(f) if m isa Product)

    # Forward-oriented reaction steps (skip inhibitor/regulator binding, pure
    # conformational isos, and product-rebinding dead-ends). For each: type ∈
    # {:bind,:chem,:release}; (ff,ft) = forward (from,to). A product release stored
    # as product-binding (`E+P→EP`, product on m_lhs) is the reverse reaction, so
    # its forward direction swaps the endpoints; an SS-dissociation release
    # (`EA→E+P`, product on m_rhs) is already forward.
    rsteps = NamedTuple[]
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        bm = bound_metabolite(s)
        local ff, ft, typ
        if bm === nothing                                   # iso step
            f, t = from_species(s), to_species(s)
            (Set(name(m) for m in bound(f)) == Set(name(m) for m in bound(t)) &&
             residual(f) == residual(t)) && continue        # pure conformational iso
            ff, ft, typ = f, t, :chem
        elseif bm isa Substrate
            ff, ft, typ = from_species(s), to_species(s), :bind
        elseif bm isa Product
            _, _, m_lhs, _ = _step_sides(s)
            rev = !isempty(m_lhs)                           # product-binding storage
            ff = rev ? to_species(s) : from_species(s)
            ft = rev ? from_species(s) : to_species(s)
            typ = :release
        else
            continue                                        # inhibitor / regulator
        end
        typ !== :chem &&
            (is_mixed_substrate_product_complex(ff) ||
             is_mixed_substrate_product_complex(ft)) && continue
        push!(rsteps, (idx = idx, ff = ff, ft = ft, typ = typ, s = s))
    end

    # Candidate cuts: (steps into rsteps, central form or `nothing` for metabolite).
    cands = Tuple{Vector{Int}, Union{Species, Nothing}}[]
    add_cand!(central, pred) = (g = [k for k in eachindex(rsteps) if pred(rsteps[k])];
                                isempty(g) || push!(cands, (g, central)))
    for S in subs_species
        add_cand!(nothing, r -> r.typ === :bind && name(bound_metabolite(r.s)) == S)
        add_cand!(nothing, r -> r.typ === :chem && S in subs_in(r.ff) && !(S in subs_in(r.ft)))
    end
    for P in prods_species
        add_cand!(nothing, r -> r.typ === :release && name(bound_metabolite(r.s)) == P)
        add_cand!(nothing, r -> r.typ === :chem && P in prods_in(r.ft) && !(P in prods_in(r.ff)))
    end
    central_forms = Set{Species}()      # iso-step endpoints
    for r in rsteps
        r.typ === :chem && (push!(central_forms, r.ff); push!(central_forms, r.ft))
    end
    for X in sort(collect(central_forms); by = x -> string(name(x)))  # sorted: precompile-stable
        add_cand!(X, r -> r.ft == X)    # produce X
        add_cand!(X, r -> r.ff == X)    # consume X
    end

    usable = [c for c in cands if all(!is_equilibrium(rsteps[k].s) for k in c[1])]
    if isempty(usable)
        allow_dead && return poly_zero()
        error(
            "rate_equation: no rapid-equilibrium-consistent reaction cut — a complete " *
            "all-RE catalytic cycle exists, so the mechanism has no finite rate.")
    end

    # Prefer a metabolite cut (central === nothing) over a central-species cut, then
    # fewest steps, then a chemistry cut, then sorted indices (precompile-stable).
    has_chem(c) = any(rsteps[k].typ === :chem for k in c[1])
    steps, central = usable[argmin(
        i -> (usable[i][2] === nothing ? 0 : 1, length(usable[i][1]),
              has_chem(usable[i]) ? 0 : 1, sort([rsteps[k].idx for k in usable[i][1]])),
        eachindex(usable))]

    # A central-species cut is complete only when its form is unique up to bound
    # regulators; a regulator-variant sibling is a parallel route this cut would
    # undercount, so raise rather than return a wrong rate.
    if central !== nothing
        xs = sort!([name(m) for m in bound(central) if m isa Substrate])
        xp = sort!([name(m) for m in bound(central) if m isa Product])
        for Y in _enumerate_species(mech)
            Y == central && continue
            sort!([name(m) for m in bound(Y) if m isa Substrate]) == xs &&
                sort!([name(m) for m in bound(Y) if m isa Product]) == xp &&
                residual(Y) == residual(central) && error(
                "rate_equation: ambiguous central-complex cut on $(name(central)) — " *
                "the regulator-variant form $(name(Y)) is a parallel route this cut " *
                "would undercount; cannot derive a unique rate.")
        end
    end

    num = poly_zero()
    for k in steps
        r = rsteps[k]; idx = r.idx
        e_lhs, e_rhs, m_lhs, m_rhs = _step_sides(r.s)
        i_form = enz_name_to_form[e_lhs]; j_form = enz_name_to_form[e_rhs]
        g1, g2 = form_to_group[i_form], form_to_group[j_form]
        rf = _ss_contrib(poly_sym(name(step_params[idx][1], mech)), m_lhs, i_form, alpha)
        rr = _ss_contrib(poly_sym(name(step_params[idx][2], mech)), m_rhs, j_form, alpha)
        canon = poly_sub(poly_mul(rf, D[g1]), poly_mul(rr, D[g2]))
        num = poly_add(num, (r.typ === :release && !isempty(m_lhs)) ? poly_neg(canon) : canon)
    end
    num
end

# ─── Expr generation from POLY ──────────────────────────────

"""
Compute the raw rate expression (bare symbols) and sorted parameter/concentration symbols.
Returns `(expr, all_params, sorted_concs)`.
"""
function _raw_rate_expr_and_symbols(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    param_syms = Set{Symbol}(_raw_param_symbols(m))
    conc_syms = Set{Symbol}(metabolites(m))
    num_expr = _poly_to_expr(num, param_syms, conc_syms)
    den_expr = _poly_to_expr(den, param_syms, conc_syms)
    expr = :(E_total * ($num_expr) / ($den_expr))
    all_params = _sorted_raw_param_symbols(M)
    return expr, all_params, metabolites(m)
end

# ─── Mode-dispatched rate_equation ────────────────────────────

"""
    rate_equation(m, concs, params, [mode])

Return the net reaction rate of mechanism `m` at the metabolite concentrations
`concs` and parameters `params` (both `NamedTuple`s). The rate equation itself
is derived from the mechanism's elementary steps by the King–Altman/Cha method:
rapid-equilibrium steps collapse to binding constants, steady-state steps are
assembled into the King–Altman determinant, and the two combine into the
quasi-steady-state flux through the whole mechanism — so one call returns the
overall turnover, not the rate of a single step.

`mode` selects how parameters enter the equation. The default [`Reduced`](@ref)
applies the Haldane/Wegscheider reduction, deriving the dependent rate constants
from `Keq` and the independent parameters; [`Full`](@ref) instead takes every
rate constant as independent.

The derivation runs once, at compile time: the body is generated as a single
arithmetic expression with no allocations, loops, or matrix operations. Use
[`rate_equation_string`](@ref) to inspect the derived equation symbolically.
"""
function rate_equation end

rate_equation(m::_AnyMechanism, concs, params) = rate_equation(m, concs, params, Reduced)

# Concrete-mechanism convenience: lift to the singleton (see the note on
# the parameters/fitted_params methods above — allocates, not the hot path).
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params,
              mode::AbstractRateEquationMode) =
    rate_equation(compile_mechanism(m), concs, params, mode)
rate_equation(m::Union{Mechanism, AllostericMechanism}, concs, params) =
    rate_equation(m, concs, params, Reduced)

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple, ::FullMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, FullMode)
end

@generated function rate_equation(
    m::M, concs::NamedTuple, params::NamedTuple,
    ::ReducedMode,
) where {M <: EnzymeMechanism}
    _build_rate_body(M, ReducedMode)
end

# ─── String Representation ────────────────────────────────────

"""
    rate_equation_string(m, [mode]) -> String

Return the symbolic rate equation for mechanism `m` as a multi-line
`String` (it returns the text — it does not print). `mode` is `Reduced`
(default) or `Full`; pass a concrete `Mechanism` / `AllostericMechanism`
or its compiled [`EnzymeMechanism`](@ref) singleton.

The string is a runnable transcript of how [`rate_equation`](@ref)
evaluates: a `(; …) = params` destructure line, a `(; …) = concs`
destructure line, then the `v = E_total * (num) / (den)` line. In
`Reduced` mode, dependent rate constants are listed first under
`# Wegscheider constraints:` and `# Haldane constraints:` headers — the
thermodynamic identities that eliminate parameters — and only the
independent set appears in the `params` destructure. In `Full` mode every
rate constant is independent, so there is no constraint section. `Full`
mode is defined for `EnzymeMechanism` only; an `AllostericEnzymeMechanism`
supports `Reduced` mode only.

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
(; K_P_E, K_S_E, k_ES_to_EP, Keq, E_total) = params
(; S, P) = concs
# Haldane constraints:
k_EP_to_ES = (1 / Keq) * K_P_E * (1 / K_S_E) * k_ES_to_EP
v = E_total * (k_ES_to_EP * S / K_S_E - k_EP_to_ES * P / K_P_E) / (1 + P / K_P_E + S / K_S_E)
```
"""
function rate_equation_string end

rate_equation_string(m::_AnyMechanism) = rate_equation_string(m, Reduced)

# Concrete-mechanism convenience: lift to the singleton.
rate_equation_string(m::Union{Mechanism, AllostericMechanism},
                     mode::AbstractRateEquationMode) =
    rate_equation_string(compile_mechanism(m), mode)
rate_equation_string(m::Union{Mechanism, AllostericMechanism}) =
    rate_equation_string(m, Reduced)

"""Build the `v = E_total * (num) / (den)` line from the raw symbolic rate polys."""
function _rate_v_line(M::Type{<:EnzymeMechanism})
    num, den = _raw_symbolic_rate_polys(M)
    m = M()
    ps = Set{Symbol}(_raw_param_symbols(m))
    cs = Set{Symbol}(metabolites(m))
    "v = E_total * ($(_expr_to_string(_poly_to_expr(num, ps, cs)))) / " *
        "($(_expr_to_string(_poly_to_expr(den, ps, cs))))"
end

function rate_equation_string(::M, ::FullMode) where {M<:EnzymeMechanism}
    mech = Mechanism(M())
    param_names = Symbol[name(p, mech) for p in _enumerate_parameters_full(mech)]
    lines = ["(; $(join((param_names..., :E_total), ", "))) = params",
             "(; $(join(metabolites(M()), ", "))) = concs"]
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end

"""
Render each dep-map entry as a constraint line and append it to `weg_lines` or
`hal_lines`, split by whether its RHS references `Keq` (Haldane) or not
(Wegscheider). Single-symbol RHSes get the substituted-into-v annotation;
multi-symbol RHSes get runtime assignment in `_build_rate_body` (no annotation).
Entries are visited in lexicographic LHS order.
"""
function _partition_constraint_lines!(weg_lines, hal_lines, dep)
    keq_set = Set([:Keq])
    for (sym, expr) in sort(collect(dep); by=p -> string(p[1]))
        is_haldane = _expr_references_any(expr, keq_set)
        suffix = expr isa Symbol ? ANNOTATION_SUBSTITUTED : ""
        push!(is_haldane ? hal_lines : weg_lines, "$sym = $(string(expr))$suffix")
    end
end

"""Append the `# Wegscheider constraints:` and `# Haldane constraints:` sections
(each skipped when empty) to `lines`."""
function _append_constraint_sections!(lines, weg_lines, hal_lines)
    isempty(weg_lines)  ||
        (push!(lines, "# Wegscheider constraints:");  append!(lines, weg_lines))
    isempty(hal_lines)  ||
        (push!(lines, "# Haldane constraints:");      append!(lines, hal_lines))
end

function rate_equation_string(::M, ::ReducedMode) where {M<:EnzymeMechanism}
    m = M()
    _, indep = _dependent_param_exprs(M)

    dep_raw, _ = _dependent_param_exprs_kernel(M, Dict{Symbol, Symbol}())
    weg_lines, hal_lines = String[], String[]
    _partition_constraint_lines!(weg_lines, hal_lines, dep_raw)

    lines = ["(; $(join((indep..., :Keq, :E_total), ", "))) = params",
             "(; $(join(metabolites(m), ", "))) = concs"]
    _append_constraint_sections!(lines, weg_lines, hal_lines)
    push!(lines, _rate_v_line(M))
    join(lines, "\n")
end

# ─── kcat Computation Helpers ──────────────────────────────────

"""
Set of Symbol names for SS rate-constant parameters (Kon, Koff, Kfor,
Krev) of `em`. For `AllostericEnzymeMechanism`, also includes the
`_T`-suffixed names of every SS rate constant that lives in the
inactive state polynomial. Routes Symbol production through the
`name(p, m)` chokepoint via Parameter-subtype dispatch. Used by
`rescale_parameter_values` to scale only SS k's without touching RE
Kd's, Keq, E_total, L, or regulatory K's.
"""
function _ss_rate_constant_names(em::EnzymeMechanism)
    mech = Mechanism(em)
    Set{Symbol}(name(p, mech) for p in _enumerate_parameters_full(mech)
                if p isa Union{Kon, Koff, Kfor, Krev})
end

function _ss_rate_constant_names(em::AllostericEnzymeMechanism)
    am = AllostericMechanism(em)
    a_names = Set{Symbol}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        rep = _group_rep(group, fes)
        is_equilibrium(rep) && continue
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        for p in _emit_cat_params_for_rep(rep, st)
            push!(a_names, name(p, am))
        end
    end
    i_names = Set{Symbol}(name(p, am) for p in _all_i_state_parameters(am)
                          if p isa Union{Kon, Koff, Kfor, Krev})
    union(a_names, i_names)
end

"""Group `num` and `den` POLYs by metabolite monomial pattern, using
`k_param_names` to classify each symbol. Returns `(num_groups,
den_groups)` where each value is a POLY of k-monomials sharing the
same met-monomial. Reverse (negative-coefficient) terms are dropped
from the numerator. Used by `_kcat_forward` to compare saturating
metabolite patterns across A/I states."""
function _kcat_groups_from_polys(num::POLY, den::POLY,
                                  k_param_names::Set{Symbol})
    function split_mono(mono::MONO)
        k_mono = MONO()
        met_mono = MONO()
        for (s, e) in mono
            if s in k_param_names || s == :Keq
                push!(k_mono, s => e)
            elseif s != :E_total
                push!(met_mono, s => e)
            end
        end
        sort!(k_mono; by=first), sort!(met_mono; by=first)
    end

    num_groups = Dict{MONO, POLY}()
    for (mono, coeff) in num
        coeff > 0 || continue
        k_part, met_part = split_mono(mono)
        p = get!(num_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    den_groups = Dict{MONO, POLY}()
    for (mono, coeff) in den
        k_part, met_part = split_mono(mono)
        p = get!(den_groups, met_part, POLY())
        p[k_part] = get(p, k_part, Rational{Int}(0)) + coeff
    end
    num_groups, den_groups
end

"""
    _kcat_forward(m::EnzymeMechanism, params) → Float64

Compute kcat (forward) analytically from the polynomial structure.
kcat is the maximum rate at saturating substrates, zero products,
and E_total=1. For mechanisms with multiple catalytic paths
(e.g., non-essential activator), returns the max over all paths.

Groups the rate-equation numerator (forward terms only, positive coefficients)
and denominator by metabolite pattern. For each matching pair, builds k-only
(num_k_expr / den_k_expr) candidates and emits `max(nk/dk, ...)`. Equilibrium
constants cancel between num and den at matching metabolite levels.

Multiple candidates arise for mechanisms with alternative catalytic pathways
(e.g., non-essential activator with/without activator bound).
"""
@generated function _kcat_forward(
    ::M, params::NamedTuple,
) where {M <: EnzymeMechanism}
    num, den = _raw_symbolic_rate_polys(M)
    k_param_names = Set{Symbol}(_raw_param_symbols(M()))
    num_groups, den_groups = _kcat_groups_from_polys(num, den, k_param_names)

    # Build kcat candidates: for each forward numerator metabolite group
    # with a matching denominator group, create (num_k_expr, den_k_expr).
    # kcat is evaluated at products = 0, so product-containing monomials are
    # outside its domain — King–Altman net-flux cross-terms like A·B·P yield
    # spurious candidates that can win the max. Keep substrate-only patterns.
    empty_set = Set{Symbol}()
    prod_syms = Set{Symbol}(products(M()))
    components = Tuple{Any, Any}[]
    for (met_key, num_k) in sort!(collect(num_groups); by=first)
        den_k = get(den_groups, met_key, nothing)
        den_k === nothing && continue
        any(first(s) in prod_syms for s in met_key) && continue
        num_expr = _poly_to_expr(num_k, empty_set, empty_set)
        den_expr = _poly_to_expr(den_k, empty_set, empty_set)
        push!(components, (num_expr, den_expr))
    end

    dep_exprs, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq)
    assignments = [Expr(:(=), sym, dep_exprs[sym])
                   for (sym, _) in sort(collect(dep_exprs); by=first)]
    candidates = [:($nk / $dk) for (nk, dk) in components]
    result = length(candidates) == 1 ?
        candidates[1] : Expr(:call, :max, candidates...)
    Expr(:block,
        _destructuring_expr(hw_params, :params),
        assignments...,
        :(return $result))
end

"""
    _kcat_forward(m::AllostericEnzymeMechanism, params) → Float64

Returns the peak achievable forward turnover: `max` over saturating
substrate patterns and regulator corners (each regulator 0 or saturating)
at products = 0, E_total = 1. Equals the numerical grid-peak forward rate.
Per active site (protomer) — `E_total` is the active-site concentration,
so this carries no `catalytic_multiplicity` factor.
"""
@generated function _kcat_forward(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    params::NamedTuple,
) where {CM,CS,RS}
    M_type = AllostericEnzymeMechanism{CM,CS,RS}
    aem = M_type()
    am  = AllostericMechanism(aem)
    CatN = catalytic_multiplicity(am)

    # Build A-state and I-state polynomials natively per conformation
    # (`_state_rate_polys`), so the saturating metabolite pattern can be matched
    # across conformations. The I-run derives on the reachable-form-pruned graph:
    # a pattern that only exists via an `:OnlyA` step drops out of I-state, and a
    # dead I-cycle yields `num_I = poly_zero()` (empty groups), so its saturating
    # contribution vanishes.
    num_A_poly, den_A_poly = _state_rate_polys(am, :A)
    num_I_poly, den_I_poly = _state_rate_polys(am, :I)
    cat_mets = Set{Symbol}(metabolites(CM()))
    # Catalytic param-name sets for the metabolite/k split. The A-set is the
    # A-state tagged column set; the I-set adds the I-polynomials' own params
    # (`:I` mirrors plus any Case-B name the native I-run introduced), which are
    # exactly the non-metabolite symbols the I-polys reference.
    a_param_names = Set(_state_all_params(_state_mechanism(am, :A),
                                          _state_step_params(am, :A)))
    i_param_names = union(a_param_names,
        setdiff(union(_poly_param_syms(num_I_poly), _poly_param_syms(den_I_poly)),
                cat_mets))
    num_A_groups, den_A_groups =
        _kcat_groups_from_polys(num_A_poly, den_A_poly, a_param_names)
    num_I_groups, den_I_groups =
        _kcat_groups_from_polys(num_I_poly, den_I_poly, i_param_names)

    # kcat = peak forward turnover at saturation: max over saturating patterns
    # (met_key) and regulator corners. Only substrate-saturating patterns are
    # valid at products=0; product-containing patterns are excluded.
    prod_syms = Set{Symbol}(name(p) for p in products(am.reaction))
    a_keys = sort!([k for k in keys(num_A_groups)
                    if haskey(den_A_groups, k) && !any(first(s) in prod_syms for s in k)])
    isempty(a_keys) &&
        error("_kcat_forward: AllostericEnzymeMechanism produced no kcat " *
              "components — saturating-substrate pattern not found in numerator")
    empty_set = Set{Symbol}()

    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # Keep inactive-state assignments unconditionally: B_I references them, and
    # deps touching an :OnlyA symbol are already zeroed in _build_dep_assignments.
    i_assignments = i_assignments_
    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq)
    i_state_dead = _i_state_num_zero(am)

    # Regulator-corner setup (independent of the saturating pattern).
    all_ligs = AllostericRegulator[]
    for site in am.regulatory_sites
        for lig in site.ligands
            lig in all_ligs || push!(all_ligs, lig)
        end
    end
    n_ligs = length(all_ligs)
    lig_idx = Dict(lig => i - 1 for (i, lig) in enumerate(all_ligs))

    kcat_exprs = Any[]
    for met_key in a_keys
        num_k_A_expr = _poly_to_expr(num_A_groups[met_key], empty_set, empty_set)
        den_k_A_expr = _poly_to_expr(den_A_groups[met_key], empty_set, empty_set)
        num_I_p = get(num_I_groups, met_key, nothing)
        den_I_p = get(den_I_groups, met_key, nothing)
        num_k_I_expr = num_I_p === nothing ? 0 :
            _poly_to_expr(num_I_p, empty_set, empty_set)
        den_k_I_expr = den_I_p === nothing ? 0 :
            _poly_to_expr(den_I_p, empty_set, empty_set)
        i_pattern_dead = den_I_p === nothing

        A_A, B_A = _mwc_power_pair(num_k_A_expr, den_k_A_expr, CatN)
        if i_pattern_dead
            A_I = 0
            B_I = 0
        else
            A_I_live, B_I = _mwc_power_pair(num_k_I_expr, den_k_I_expr, CatN)
            A_I = i_state_dead ? 0 : A_I_live
        end

        if isempty(RS)
            push!(kcat_exprs,
                  :($(_mwc_combine(A_A, A_I)) / $(_mwc_combine(B_A, B_I))))
        else
            for mask in 0:(2^n_ligs - 1)
                W_A_factors = Any[]
                W_I_factors = Any[]
                for (site_idx, site) in enumerate(am.regulatory_sites)
                    n_reg = site.multiplicity
                    sat_terms_A = Any[]
                    sat_terms_I = Any[]
                    for (lig, tag) in zip(site.ligands, site.allo_states)
                        if (mask >> lig_idx[lig]) & 1 == 1
                            if tag !== :OnlyI
                                K_A_sym = name(Kreg(site, lig, :A), am)
                                push!(sat_terms_A, :(inv($K_A_sym)))
                            end
                            if tag !== :OnlyA
                                # `:EqualAI` ligands share the A-state symbol; the
                                # body emits an `:EqualAI` ligand's I-state slot
                                # via the A-state name (no `_T` rename).
                                K_I_state = tag === :EqualAI ? :A : :I
                                K_I_sym = name(Kreg(site, lig, K_I_state), am)
                                push!(sat_terms_I, :(inv($K_I_sym)))
                            end
                        end
                    end
                    if !isempty(sat_terms_A)
                        q_A = length(sat_terms_A) == 1 ?
                            sat_terms_A[1] : _nest_binary(:+, sat_terms_A)
                        push!(W_A_factors, _power_expr(q_A, n_reg))
                    end
                    if !isempty(sat_terms_I)
                        q_I = length(sat_terms_I) == 1 ?
                            sat_terms_I[1] : _nest_binary(:+, sat_terms_I)
                        push!(W_I_factors, _power_expr(q_I, n_reg))
                    end
                end
                if isempty(W_A_factors) && isempty(W_I_factors)
                    kcat_expr =
                        :($(_mwc_combine(A_A, A_I)) / $(_mwc_combine(B_A, B_I)))
                else
                    W_A = isempty(W_A_factors) ? 1 :
                        length(W_A_factors) == 1 ? W_A_factors[1] :
                        _nest_binary(:*, W_A_factors)
                    W_I = isempty(W_I_factors) ? 1 :
                        length(W_I_factors) == 1 ? W_I_factors[1] :
                        _nest_binary(:*, W_I_factors)
                    kcat_expr = :(($(A_A) * $(W_A) + L * $(A_I) * $(W_I)) /
                        ($(B_A) * $(W_A) + L * $(B_I) * $(W_I)))
                end
                push!(kcat_exprs, kcat_expr)
            end
        end
    end

    result = length(kcat_exprs) == 1 ? kcat_exprs[1] :
        Expr(:call, :max, kcat_exprs...)
    return Expr(:block,
        _destructuring_expr(hw_params, :params),
        a_assignments...,
        i_assignments...,
        :(return $result))
end

# ─── Public API: rescale_parameter_values ──────────────────────────

"""
    rescale_parameter_values(m, params::NamedTuple; scale_k_to_kcat=1.0)

Rescale SS rate constants so that `_kcat_forward(m, result) ≈ scale_k_to_kcat`.
Non-SS parameters (K's, Keq, E_total, L, regulatory K's) are unchanged.
"""
function rescale_parameter_values(
    m::_AnyMechanism, params::NamedTuple; scale_k_to_kcat=1.0,
)
    kcat_current = _kcat_forward(m, params)
    scale = scale_k_to_kcat / kcat_current
    ss_names = _ss_rate_constant_names(m)
    NamedTuple{keys(params)}(Tuple(
        k in ss_names ? v * scale : v
        for (k, v) in zip(keys(params), values(params))
    ))
end

# ═══════════════════════════════════════════════════════════════════
# AllostericEnzymeMechanism rate equations (MWC)
#
# Per-active-site (per-protomer) normalization: `E_total` is the active-
# site concentration, so the leading `cat_n` multiplicity coefficient is
# absorbed into E_total and does not appear in the rate. The
# `Q_cat_c^(cat_n - 1)` / `Q_cat_c^cat_n` binding-statistics powers stay.
#
# MWC rate formula (per conformation c, summed over conformations).
# Let cat_n = catalytic_multiplicity(m):
#   num = sum_c( L_c * N_cat_c * Q_cat_c^(cat_n - 1)
#             * prod(Q_reg_i_c^n_reg_i for all regulatory sites i) )
#   den = sum_c( L_c * Q_cat_c^cat_n * prod(Q_reg_i_c^n_reg_i for all regulatory sites i) )
#   v = E_total * num / den
#
# Regulatory sites contribute to BOTH numerator and denominator at their
# multiplicity, regardless of whether n_reg_i matches cat_n.
# ═══════════════════════════════════════════════════════════════════

"""
The native I-state catalytic numerator polynomial is empty (zero): the
reachable-form-pruned I-graph of a broken cycle has no steady-state cut, so
`_compute_numerator(allow_dead=true)` returns `poly_zero()` natively — no forced
zero. It decides, at all three consumer sites (`_allosteric_num_den_exprs`,
`_i_state_referenced_syms`, `_kcat_forward`), whether the `L·num_I` term is
emitted, whether `num_I`'s symbols are defined, and whether `kcat` carries the
I-state term. A live redundant-path `:OnlyA` mechanism (num_I ≠ 0) is thereby
handled consistently.
"""
_i_state_num_zero(am::AllostericMechanism) =
    isempty(first(_state_rate_polys(am, :I)))

"""
Names of enzyme forms in the connected component of a free (empty-bound) form
over ALL steps of `groups` (rapid-equilibrium and steady-state alike). Seeding
from every empty-bound form covers a ping-pong covalent intermediate, which
carries no bound metabolite and so is its own free root.
"""
function _reachable_from_free(groups)
    forms = Species[]
    for grp in groups, s in grp
        from_species(s) in forms || push!(forms, from_species(s))
        to_species(s)   in forms || push!(forms, to_species(s))
    end
    reach = Set{Symbol}(name(f) for f in forms if isempty(bound(f)))
    changed = true
    while changed
        changed = false
        for grp in groups, s in grp
            fn, tn = name(from_species(s)), name(to_species(s))
            if fn in reach && tn ∉ reach
                push!(reach, tn); changed = true
            elseif tn in reach && fn ∉ reach
                push!(reach, fn); changed = true
            end
        end
    end
    reach
end

"""
The `AllostericMechanism` for `am` in conformational `state`: `am` itself for
`:A`; for `:I`, a fresh `AllostericMechanism` with `:OnlyA` catalytic groups
dropped AND every enzyme form disconnected from free E by that drop pruned at
the step level. After removing the `:OnlyA` groups, a form is kept iff it lies
in the connected component of a free (empty-bound) enzyme form over ALL
remaining steps (rapid-equilibrium and steady-state alike); a step is kept iff
both its endpoints are kept, so a kinetic group with all its steps dropped
disappears. Forms whose only route back to free E ran through an `:OnlyA` group
become disconnected and drop out; forms still reachable through the surviving
steps — including a substrate complex repopulated by reverse catalysis — are
retained. A mechanism with no `:OnlyA` catalytic group keeps its whole graph
and re-derives its full native I-state.

Routing both the derivation mechanism and the step_params through this one
struct keeps steps and allo-state tags aligned: the `AllostericMechanism`
constructor canonicalizes `cat_steps` and applies the SAME permutation to
`cat_allo_states`, so pruning-induced reorder/iso-flips can't desync the two.
"""
function _state_allo_mechanism(am::AllostericMechanism, state::Symbol)
    state === :I || return am
    keepG = [g for g in eachindex(steps(am)) if cat_allo_state(am, g) !== :OnlyA]
    groups = steps(am)[keepG]
    states = cat_allo_states(am)[keepG]
    all_forms = Set{Symbol}()
    for grp in groups, s in grp
        push!(all_forms, name(from_species(s)), name(to_species(s)))
    end
    stranded = setdiff(all_forms, _reachable_from_free(groups))
    kg = Vector{Step}[]
    ks = Symbol[]
    for (grp, st) in zip(groups, states)
        kept = [s for s in grp
                if name(from_species(s)) ∉ stranded &&
                   name(to_species(s)) ∉ stranded]
        isempty(kept) || (push!(kg, kept); push!(ks, st))
    end
    AllostericMechanism(reaction(am), kg, ks,
                        catalytic_multiplicity(am), regulatory_sites(am))
end

"""
State-tagged `step_params` for `am`'s catalytic mechanism in conformational
`state` (`:A` or `:I`). Same shape as `_step_parameters(::Mechanism)` — a
per-flat-step vector of `Parameter`s — but each catalytic group's `Parameter`s
carry the group's state tag: a `:NonequalAI`/`:OnlyA` group is tagged with
`state`; an `:EqualAI` group is tagged `:EqualAI` (so `name(p, am)` renders the
shared bare Symbol in both states). For `state == :I`, `:OnlyA` groups are
already pruned from `_state_allo_mechanism(am, :I)`, matching the broken-cycle
graph from `_state_mechanism(am, :I)`.

Walks `_state_allo_mechanism`'s already-canonical, aligned steps/states so the
per-flat-step order matches `_flat_steps(_state_mechanism(am, state))`. Each
`Parameter` is anchored on its own step (not the rep), exactly as
`_step_parameters` does, so arity follows the step's RE/SS type while
`name(p, am)` collapses to the rep's structural Symbol via the chokepoint.
"""
function _state_step_params(am::AllostericMechanism, state::Symbol)
    sam = _state_allo_mechanism(am, state)
    out = Vector{Vector{Parameter}}()
    for (g, group) in enumerate(steps(sam))
        tag = cat_allo_state(sam, g) === :EqualAI ? :EqualAI : state
        for s in group
            push!(out, is_equilibrium(s) ?
                Parameter[is_binding(s) ? Kd(s, tag) : Kiso(s, tag)] :
                Parameter[is_binding(s) ? Kon(s, tag)  : Kfor(s, tag),
                          is_binding(s) ? Koff(s, tag) : Krev(s, tag)])
        end
    end
    out
end

"""
The catalytic `Mechanism` for `am` in conformational `state`. `:A` keeps the
full catalytic graph; `:I` keeps only the reachable-form subgraph (`:OnlyA`
groups and the forms they disconnect from free E pruned, via
`_state_allo_mechanism`) so King–Altman re-derives the broken-cycle I-state law
natively.
"""
function _state_mechanism(am::AllostericMechanism, state::Symbol)
    sam = _state_allo_mechanism(am, state)
    Mechanism(reaction(sam), steps(sam))
end

"""
Derive `(num_poly, den_poly)` for `am`'s catalytic mechanism in conformational
`state`, natively in that state's parameter names. Runs the shared King–Altman
engine on the state-tagged `step_params` and state graph, so no post-hoc rename
is needed (`:EqualAI` groups render the shared bare Symbol automatically). For
`:I`, applies the one-rule Case-B naming (a shared `:EqualAI` dependent whose
Haldane RHS references a `:NonequalAI` symbol takes its distinct I-name) so the
polynomials reference the same I-symbols the dep-assignment preamble defines.
"""
function _state_rate_polys(am::AllostericMechanism, state::Symbol)
    cm = _state_mechanism(am, state)
    sp = _state_step_params(am, state)
    @assert length(sp) == length(_flat_steps(cm)) "state step_params/steps misaligned"
    subs_syms = Symbol[name(s) for s in substrates(reaction(am))]
    prods_syms = Symbol[name(p) for p in products(reaction(am))]
    num, den = _raw_symbolic_rate_polys(cm, sp,
                                        _state_wegscheider_rename_map(am, state),
                                        subs_syms, prods_syms;
                                        allow_dead = state === :I)
    state === :I || return num, den
    renames = _state_i_case_b_renames(am)
    _rename_symbols(num, renames), _rename_symbols(den, renames)
end

"""
Tagged catalytic parameter symbols (the kernel's column set) for state graph
`cm` under the state-tagged `sp`, in group order — the state-tagged analog of
`_raw_param_symbols`. Distinct `name(p, cm)` (rep names, no Wegscheider rename);
the kernel applies `_state_wegscheider_rename_map` on top of this column set when
folding single-symbol RE binding-K ties.
"""
function _state_all_params(cm::Mechanism, sp)
    out = Symbol[]
    seen = Set{Symbol}()
    for group in sp, p in group
        s = name(p, cm)
        s in seen || (push!(seen, s); push!(out, s))
    end
    out
end

"""
State-tagged Wegscheider rename map for `am`'s catalytic sub-mechanism in
conformational `state` — the state-aware analog of `_build_wegscheider_rename_map`
(which runs the kernel with `:None` step_params and so cannot see the `:A`/`:I`
tags). Discovers single-symbol RE binding-K Wegscheider ties (`K_a = K_b`, both
binding K's) under the state-tagged `step_params`/`all_params` and folds each
absorbed symbol into its target. Empty for all current specs (catalysis is
steady-state, so no fully-RE catalytic box), but a fully-RE catalytic core would
now collapse its tie natively — the same way the non-allosteric path does.
"""
function _state_wegscheider_rename_map(am::AllostericMechanism, state::Symbol)
    cm = _state_mechanism(am, state)
    sp = _state_step_params(am, state)
    rename = Dict{Symbol, Symbol}()
    # binding-K set: value-context rep name of each RE binding step.
    binding_set = Set{Symbol}()
    for (idx, (s, _)) in enumerate(_flat_steps(cm))
        is_equilibrium(s) && is_binding(s) || continue
        push!(binding_set, name(sp[idx][1], cm))
    end
    # Single-symbol Wegscheider RE ties between two binding K's.
    dep_raw, _ = _dependent_param_exprs_kernel(cm, rename;
                                               step_params = sp,
                                               all_params = _state_all_params(cm, sp))
    for (lhs, rhs) in dep_raw
        rhs isa Symbol || continue
        lhs in binding_set && rhs in binding_set || continue
        target = get(rename, rhs, rhs)
        rename[lhs] = target
        for k in collect(keys(rename))
            rename[k] == lhs && (rename[k] = target)
        end
    end
    rename
end

"""
Native per-state Haldane/Wegscheider dependent-parameter expressions from the
shared kernel, in that state's parameter names (NO Case-B rename yet). Runs the
kernel under the state-tagged `step_params`/`all_params` with
`_state_wegscheider_rename_map` so a fully-RE catalytic binding-K Wegscheider tie
is collapsed natively (the same absorption the non-allosteric
`_dependent_param_exprs(CM)` wrapper applies). Absorbed symbols are filtered from
`indep` (they no longer appear in `v` and must not be fittable dummies). Empty
rename for all current specs (catalysis is steady-state), so byte-identical
there.
"""
function _state_raw_dependent_exprs(am::AllostericMechanism, state::Symbol)
    cm = _state_mechanism(am, state)
    sp = _state_step_params(am, state)
    @assert length(sp) == length(_flat_steps(cm)) "state step_params/steps misaligned"
    rename = _state_wegscheider_rename_map(am, state)
    dep, indep = _dependent_param_exprs_kernel(cm, rename;
                                               step_params = sp,
                                               all_params = _state_all_params(cm, sp))
    indep = Tuple(p for p in indep if get(rename, p, p) == p)
    dep, indep
end

"""
I-state catalytic parameter Symbols of `:NonequalAI` groups — the names that
genuinely differ between conformations (`K_I_…`/`k_I_…`), emitted through the
`name(p, am)` chokepoint. Marks which symbols a Case-B dependent's RHS must
reference to earn its own I-name.
"""
function _i_nonequalai_syms(am::AllostericMechanism)
    out = Set{Symbol}()
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :NonequalAI || continue
        for p in _emit_cat_params_for_rep(_group_rep(group, fes), :I)
            push!(out, name(p, am))
        end
    end
    out
end

"""
Case-B rename map for the native I-run (Symbol → Symbol). A dependent whose key
is a bare/`:EqualAI` symbol but whose derived RHS references an I-tagged
`:NonequalAI` symbol has a genuinely different I-value and needs a distinct
I-name (e.g. PK `k_EATPPyruvate_to_EADPPEP → k_I_EATPPyruvate_to_EADPPEP`). Deps
already carrying the I-tag (a `:NonequalAI` dep) are left alone; since Gaussian
elimination expresses each dependent purely in terms of independents, this is
one non-transitive pass.
"""
function _case_b_rename_map(dep, am::AllostericMechanism)
    i_nonequalai = _i_nonequalai_syms(am)
    renames = Dict{Symbol, Symbol}()
    isempty(i_nonequalai) && return renames
    for (k, v) in dep
        k in i_nonequalai && continue
        _expr_references_any(v, i_nonequalai) || continue
        renames[k] = _dep_inactive_name(am, k)
    end
    renames
end

"""Case-B I-rename map for `am`, derived from the native I-run deps."""
_state_i_case_b_renames(am::AllostericMechanism) =
    _case_b_rename_map(first(_state_raw_dependent_exprs(am, :I)), am)

"""
Native per-state dependent-parameter assignments `(dep_exprs, indep)` in that
state's parameter names, via the shared kernel on the state graph. For `:I`,
applies the one-rule Case-B naming so a shared `:EqualAI` dependent whose value
differs between states gets its distinct I-name (Spec §4/§4a).
"""
function _state_dependent_exprs(am::AllostericMechanism, state::Symbol)
    dep, indep = _state_raw_dependent_exprs(am, state)
    state === :I || return dep, indep
    renames = _case_b_rename_map(dep, am)
    isempty(renames) && return dep, indep
    renamed = Dict{Symbol, Union{Symbol, Expr}}()
    for (k, v) in dep
        renamed[get(renames, k, k)] = v
    end
    renamed, indep
end

"""
All I-state `Parameter`s the rate-equation body emits as constraint LHSes
— `Parameter` form. Catalytic groups: every non-`:OnlyA` group
contributes I-state Parameter(s) for its rep step (`Kd`/`Kiso`/`Kon`+
`Koff`/`Kfor`+`Krev`). Regulator sites: every non-`:OnlyA` ligand
contributes an I-state `Kreg`. Synthesized-dep I-mirrors (deps whose RHS
references a `:NonequalAI` symbol) belong to dep-parameter machinery and
are emitted Symbol-level by the dep-assignment builder.
"""
function _all_i_state_parameters(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        cat_allo_state(am, g) === :OnlyA && continue
        append!(out, _emit_cat_params_for_rep(_group_rep(group, fes), :I))
    end
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyA && continue
            push!(out, Kreg(site, lig, :I))
        end
    end
    out
end

"""
Enumerate every raw rate-constant `Parameter` for an `AllostericMechanism`
(the complete Full-mode symbol set). Order is:

1. Catalytic A-state Parameter per kinetic group (every group). `:OnlyA`
   and `:NonequalAI` use `state = :A`; `:EqualAI` uses `state = :EqualAI`
   because the symbol is shared with the I-state branch (the chokepoint
   `name(p, m)` renders both to the same `Symbol`).
2. Catalytic I-state mirrors via `_all_i_state_parameters` (skips
   `:OnlyA` groups).
3. Reg-site `Kreg(site, lig, :I)` and `Kreg(site, lig, :A)` per ligand
   (the I-state set skips `:OnlyA` ligands; the A-state set skips
   `:OnlyI` ligands).
4. `Lallo()` for the MWC coupling `L`.

Synthesized-dep I-symbols (derived deps whose RHS references a
`:NonequalAI` symbol) are NOT included — they have no Parameter
representation and are emitted Symbol-level by the caller.

Non-appearing names (e.g., a catalytic I-state mirror that's elided in
a `t_state_dead` mechanism) are harmless to include: unused names are
inert — nothing consumes a name that does not appear in the
rate-equation Exprs. This enumeration intentionally over-emits.
"""
function _enumerate_parameters_full_allosteric(am::AllostericMechanism)
    out = Parameter[]
    fes = _free_enz_set(am)
    for (g, group) in enumerate(steps(am))
        rep = _group_rep(group, fes)
        st = cat_allo_state(am, g) === :EqualAI ? :EqualAI : :A
        append!(out, _emit_cat_params_for_rep(rep, st))
    end
    append!(out, _all_i_state_parameters(am))
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyI && continue
            push!(out, Kreg(site, lig, :A))
        end
    end
    push!(out, Lallo())
    out
end

# ─── Dependent parameter expressions ─────────────────────────────

# Distinct inactive-state name for a *dependent* parameter being promoted to
# per-state (Case B: an `:EqualAI` dep whose Haldane/Wegscheider RHS references
# a `:NonequalAI` symbol). For a `:NonequalAI`/`:A` dep, `_flip_to_inactive`
# already yields a distinct `:I` name; for an `:EqualAI` dep it is a no-op, so
# fall back to the forced `:I` variant to avoid a self-map.
function _dep_inactive_name(am, k::Symbol)
    p = _param_for_symbol(am, k)
    nm = name(_flip_to_inactive(p), am)
    nm == k ? name(_force_inactive(p), am) : nm
end

"""
I-state parameter Symbols actually referenced by the retained rate-equation
polynomials, sourced from the NATIVE per-state I-polynomials: `den_i` always
(`Q_I` is kept as enzyme mass), plus `num_i` when the I-state cycle is live (a
dead cycle's native `num_i` is `poly_zero()`, contributing no symbols). This is
the single source of truth for which I-state names get defined; the `isempty`
gate is the `_i_state_num_zero` native check.
"""
function _i_state_referenced_syms(am::AllostericMechanism)
    num_i, den_i = _state_rate_polys(am, :I)
    S = _poly_param_syms(den_i)
    isempty(num_i) || union!(S, _poly_param_syms(num_i))
    S
end

"""
    _dependent_param_exprs(M::Type{<:AllostericEnzymeMechanism})

Return `(dep_exprs, indep_params)` for an AllostericEnzymeMechanism from the
NATIVE per-state derivations. A-state entries come from
`_state_dependent_exprs(am, :A)`; I-state entries from
`_state_dependent_exprs(am, :I)`, kept only for names a retained I-polynomial
references (`_i_state_referenced_syms`, `S_I`). `:EqualAI` reg mirrors
(`K_I_reg = K_A_reg`) are the only reg entries added to the dep map. Reg-site
`Kreg` names and `L` complete the independent set.

`indep` order is content-canonical by concatenated component: native A-indep
(order matches the non-allosteric catalytic derivation), then the genuinely-I
independents referenced by the polynomials, then A-state reg Kregs (canonical
`regulatory_sites` order), then I-state reg Kregs, then `L`.
"""
function _dependent_param_exprs(
    ::Type{AllostericEnzymeMechanism{CM,CS,RS}},
) where {CM,CS,RS}
    am = AllostericMechanism(AllostericEnzymeMechanism{CM,CS,RS}())
    dep_A, indep_A = _state_dependent_exprs(am, :A)
    dep_I, indep_I = _state_dependent_exprs(am, :I)
    S_I = _i_state_referenced_syms(am)

    dep = Dict{Symbol, Union{Symbol, Expr}}(dep_A)
    # I-state deps whose LHS a retained polynomial references (`Q_I` always,
    # plus `N_I` when the I-cycle is live). Covers `:NonequalAI` Case-A deps
    # and `:EqualAI` Case-B synthesized I-deps (the native I-run named them).
    for (k, v) in dep_I
        k in S_I && (dep[k] = v)
    end

    # I-state independents that are (a) genuinely distinct from the A-state
    # symbol — an `:EqualAI` group shares its bare symbol with A and is already
    # in `indep_A` — and (b) referenced by a retained I-polynomial.
    a_set = Set(indep_A)
    indep_I_list = Symbol[p for p in indep_I if p ∉ a_set && p in S_I]

    # Reg-site Parameters via `Kreg` structs + the `name(::Kreg, am)`
    # chokepoint. `:EqualAI` reg ligands share their value across states, so the
    # I-name is a dependent mirror of the A-name; `:NonequalAI`/`:OnlyI` carry a
    # distinct independent I-name.
    reg_params_a = Symbol[]
    reg_params_i_indep = Symbol[]
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :OnlyI || push!(reg_params_a, name(Kreg(site, lig, :A), am))
            if tag === :EqualAI
                dep[name(Kreg(site, lig, :I), am)] = name(Kreg(site, lig, :A), am)
            elseif tag === :NonequalAI || tag === :OnlyI
                push!(reg_params_i_indep, name(Kreg(site, lig, :I), am))
            end
        end
    end

    merged_indep = (indep_A..., indep_I_list...,
                    reg_params_a..., reg_params_i_indep..., :L)
    return dep, merged_indep
end

# `parameters` and `fitted_params` for `AllostericEnzymeMechanism`
# dispatch on explicit per-type methods at the top of this file.

# ─── Rate body building helpers ───────────────────────────────────

"""Build the regulatory site partition function expression: 1 + lig/K_lig_reg_i + ...
Skips ligands absent from the given conformation (`:OnlyA` in I-state, `:OnlyI`
in A-state). Uses the A-state K symbol when the ligand tag is `:EqualAI`.
Renders K-names via the `name(::Kreg, am)` chokepoint."""
function _reg_site_expr(am::AllostericMechanism, site_idx::Int, inactive::Bool)
    site = regulatory_sites(am)[site_idx]
    terms = Any[1]
    for (lig, tag) in zip(ligands(site), allo_states(site))
        if inactive
            tag === :OnlyA && continue
        else
            tag === :OnlyI && continue
        end
        # `:EqualAI` ligands share the A-state symbol in both conformations;
        # `:NonequalAI` / `:OnlyI` ligands carry a distinct I-state K name.
        state = (inactive && tag in (:NonequalAI, :OnlyI)) ? :I : :A
        K_sym = name(Kreg(site, lig, state), am)
        push!(terms, :($(name(lig)) / $K_sym))
    end
    _nest_binary(:+, terms)
end

"""Raise an expression to an integer power (returns 1 for n=0, expr for n=1)."""
function _power_expr(expr, n::Int)
    n == 0 && return 1
    n == 1 && return expr
    :(($expr)^$n)
end

"""MWC active + L·inactive state-combine `A + L * B`. Shared by `_kcat_forward`
(numerator/denominator halves of the state ratio) and `_allosteric_num_den_exprs`
(the retained num/den sums)."""
_mwc_combine(a, b) = :($a + L * $b)

"""MWC binding-statistics power-pair `(X * Y^(n-1), Y^n)` for a saturating pattern:
the numerator carries one fewer denominator power than the denominator, where
`n = catalytic_multiplicity`. Used by `_kcat_forward` per conformation."""
_mwc_power_pair(x, y, n) =
    (n == 1 ? x : :($x * $y^$(n - 1)), :($y^$n))

"""
Build active-state and inactive-state dep-param assignment Exprs from the NATIVE
per-state derivations. Returns `(a_assignments::Vector{Expr},
i_assignments::Vector{Expr})`. Shared by `_build_allosteric_rate_body` and
`rate_equation_string`.

A-assignments are the native A-state deps; I-assignments are the `:EqualAI` reg
mirrors (`K_I_reg = K_A_reg`) plus the native I-state catalytic deps FILTERED to
`_i_state_referenced_syms` (so no assignment is emitted for a symbol the retained
I-polynomials never reference). `:EqualAI` catalytic params share their bare
symbol across states and so need no I-mirror. All Symbols route through the
`name(p, am)` chokepoint (native derivation + `Kreg`).
"""
function _build_dep_assignments(
    M_type::Type{<:AllostericEnzymeMechanism},
)
    am = AllostericMechanism(M_type())

    dep_A, _ = _state_dependent_exprs(am, :A)
    dep_I, _ = _state_dependent_exprs(am, :I)
    # An I-state catalytic dep is emitted iff its LHS is referenced by a retained
    # I-polynomial (`Q_I` always, plus `N_I` when the I-cycle is live).
    i_names_set = _i_state_referenced_syms(am)

    a_assignments = Expr[Expr(:(=), sym, rhs)
                         for (sym, rhs) in sort(collect(dep_A); by=first)]

    i_assignments = Expr[]

    # `:EqualAI` reg params: K_I_reg = K_A_reg (before any dep that reads them).
    # Routes through the `Kreg` chokepoint.
    for site in regulatory_sites(am)
        for (lig, tag) in zip(ligands(site), allo_states(site))
            tag === :EqualAI || continue
            push!(i_assignments,
                  Expr(:(=), name(Kreg(site, lig, :I), am),
                             name(Kreg(site, lig, :A), am)))
        end
    end

    # Native I-state catalytic deps, S_I-filtered. Covers `:NonequalAI` Case-A
    # deps and `:EqualAI` Case-B synthesized I-deps the native I-run named.
    for (sym, rhs) in sort(collect(dep_I); by=first)
        sym in i_names_set && push!(i_assignments, Expr(:(=), sym, rhs))
    end

    return a_assignments, i_assignments
end

"""
Assemble the MWC numerator and denominator Exprs.
Returns `(full_num, full_den)`. Per-active-site normalization: the
numerator carries no leading `catalytic_multiplicity` factor; only the
`Q_cat^(CatN-1)` / `Q_cat^CatN` binding-statistics powers remain.
"""
function _allosteric_num_den_exprs(M_type::Type{<:AllostericEnzymeMechanism})
    m = M_type()
    am = AllostericMechanism(m)
    CM = typeof(catalytic_mechanism(m))
    CatN = catalytic_multiplicity(m)
    RS = regulatory_sites(am)

    num_A_poly, den_A_poly = _state_rate_polys(am, :A)
    # A-state catalytic param symbols (the tagged column set) drive `_poly_to_expr`'s
    # param/metabolite ordering split; the I-poly's `:I` symbols sort as non-params.
    cat_params = Set(_state_all_params(_state_mechanism(am, :A),
                                       _state_step_params(am, :A)))
    cat_mets = Set{Symbol}(metabolites(CM()))

    N_A = _poly_to_expr(num_A_poly, cat_params, cat_mets)
    Q_A = _poly_to_expr(den_A_poly, cat_params, cat_mets)

    # I-state catalytic Exprs, always re-derived natively on the reachable-form
    # subgraph (`_state_allo_mechanism(am, :I)` drops `:OnlyA` groups and every
    # form they disconnect from free E). Reachable-subgraph King–Altman gives the
    # same binding partition monomial-zeroing produced, and for a dead cycle the
    # pruned graph has no SS cut so `_compute_numerator(allow_dead=true)` returns
    # 0 natively — no forced zero needed.
    num_i_poly, den_i_poly = _state_rate_polys(am, :I)
    N_I = _poly_to_expr(num_i_poly, cat_params, cat_mets)
    Q_I = _poly_to_expr(den_i_poly, cat_params, cat_mets)

    reg_Q_A = Any[_reg_site_expr(am, i, false) for i in eachindex(RS)]
    reg_Q_I = Any[_reg_site_expr(am, i, true) for i in eachindex(RS)]

    # Numerator: N × Q_cat^(CatN-1) × all reg-site factors at multiplicity.
    function make_num_term(N, Q, reg_Qs)
        factors = Any[N]
        CatN > 1 && push!(factors, _power_expr(Q, CatN - 1))
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i], multiplicity(RS[i])))
        end
        _nest_binary(:*, factors)
    end

    # Denominator: Q_cat^CatN × all reg-site factors at multiplicity.
    function make_den_term(Q, reg_Qs)
        factors = Any[_power_expr(Q, CatN)]
        for i in eachindex(RS)
            push!(factors, _power_expr(reg_Qs[i], multiplicity(RS[i])))
        end
        _nest_binary(:*, factors)
    end

    num_A = make_num_term(N_A, Q_A, reg_Q_A)
    den_A = make_den_term(Q_A, reg_Q_A)
    den_I = make_den_term(Q_I, reg_Q_I)
    full_den = _mwc_combine(den_A, den_I)

    if isempty(num_i_poly)
        # Native I-state numerator is zero (`_i_state_num_zero`): the I-state
        # cycle is dead, so drop the L*num_I term entirely (skip dead numerator
        # branch). Q_I still contributes to denominator as enzyme mass.
        num_A, full_den
    else
        num_I = make_num_term(N_I, Q_I, reg_Q_I)
        _mwc_combine(num_A, num_I), full_den
    end
end

"""Build the MWC rate equation body as an Expr block."""
function _build_allosteric_rate_body(M_type::Type{<:AllostericEnzymeMechanism})
    full_num, full_den = _allosteric_num_den_exprs(M_type)
    rate_expr = :(E_total * ($full_num) / ($full_den))

    a_assignments, i_assignments_ = _build_dep_assignments(M_type)
    # Keep inactive-state assignments unconditionally: the retained Q_I
    # (`L * den_I`) references them. Deps whose RHS touches an :OnlyA symbol
    # are already zeroed in `_build_dep_assignments`, so nothing is undefined.
    i_assignments = i_assignments_

    _, indep = _dependent_param_exprs(M_type)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(M_type())

    Expr(:block,
        _destructuring_expr(hw_params, :params),
        _destructuring_expr(mets, :concs),
        a_assignments...,
        i_assignments...,
        :(return $rate_expr))
end

# ─── Rate equation dispatch ───────────────────────────────────────

@generated function rate_equation(
    ::AllostericEnzymeMechanism{CM,CS,RS},
    concs::NamedTuple, params::NamedTuple, ::ReducedMode,
) where {CM,CS,RS}
    _build_allosteric_rate_body(AllostericEnzymeMechanism{CM,CS,RS})
end

# ─── String representation ────────────────────────────────────────

function rate_equation_string(
    ::AllostericEnzymeMechanism{CM,CS,RS}, ::ReducedMode,
) where {CM,CS,RS}
    M = AllostericEnzymeMechanism{CM,CS,RS}
    m = M()
    am = AllostericMechanism(m)
    _, indep = _dependent_param_exprs(M)
    hw_params = (indep..., :Keq, :E_total)
    mets = metabolites(m)

    # Active-state Wegscheider/Haldane from the native A-state derivation, in
    # `:A`-state names directly (no post-hoc rename).
    dep_A, _ = _state_dependent_exprs(am, :A)
    weg_lines, hal_lines = String[], String[]
    _partition_constraint_lines!(weg_lines, hal_lines, dep_A)

    # Inactive-state assignments — partitioned by Keq-reference predicate so
    # they fold into the same two sections as active-state lines. Allosteric
    # mechanisms thus use the same three-section structure as
    # non-allosteric.
    keq_set = Set([:Keq])
    _, i_assignments_ = _build_dep_assignments(M)
    i_assignments = i_assignments_
    for a in i_assignments
        sym = a.args[1]
        expr = a.args[2]
        is_haldane = _expr_references_any(expr, keq_set)
        line = "$sym = $(_expr_to_string(expr))"
        push!(is_haldane ? hal_lines : weg_lines, line)
    end

    # Sort each section lexicographically — load-bearing for eq_hash
    # dedup of allosteric Source-C clusters, since inactive-state lines are
    # appended in iteration order rather than the lexicographic active-state
    # order.
    sort!(weg_lines)
    sort!(hal_lines)

    full_num, full_den = _allosteric_num_den_exprs(M)
    v_line = "v = E_total * ($(_expr_to_string(full_num))) / ($(_expr_to_string(full_den)))"

    lines = ["(; $(join(hw_params, ", "))) = params",
             "(; $(join(mets, ", "))) = concs"]
    _append_constraint_sections!(lines, weg_lines, hal_lines)
    push!(lines, v_line)
    join(lines, "\n")
end

"""
Canonical kinetic-group partition: merge kinetic groups whose binding-K
representatives are single-symbol Wegscheider-tied (the relation
`_build_wegscheider_rename_map` finds), so split and merged encodings of the
same rate-equivalent graph collapse to one partition. Returns `mech.steps`
unchanged when nothing is tied.
"""
function _merge_tied_kinetic_groups(mech::Mechanism)
    rename = _build_wegscheider_rename_map(mech)
    isempty(rename) && return mech.steps
    groups = mech.steps
    step_params = _step_parameters(mech)
    grp_of_flat = Int[]
    for (gi, g) in enumerate(groups), _ in g
        push!(grp_of_flat, gi)
    end
    # Canonical representative binding-K name per group (nothing if no binding step).
    rep = Vector{Union{Symbol, Nothing}}(nothing, length(groups))
    for (idx, (s, _)) in enumerate(_flat_steps(mech))
        is_equilibrium(s) && is_binding(s) || continue
        k = name(step_params[idx][1], mech)
        rep[grp_of_flat[idx]] = get(rename, k, k)
    end
    byrep = Dict{Symbol, Vector{Int}}()
    for (gi, r) in enumerate(rep)
        r === nothing && continue
        push!(get!(byrep, r, Int[]), gi)
    end
    any(length(v) > 1 for v in values(byrep)) || return mech.steps
    merged = Vector{Vector{Step}}()
    done = falses(length(groups))
    for gi in eachindex(groups)
        done[gi] && continue
        r = rep[gi]
        if r !== nothing && length(byrep[r]) > 1
            gis = byrep[r]
            push!(merged, Step[s for j in gis for s in groups[j]])
            for j in gis
                done[j] = true
            end
        else
            push!(merged, copy(groups[gi]))
            done[gi] = true
        end
    end
    merged
end
