# ABOUTME: Runs the slice's three enumeration moves on every parent set the package run saved,
# ABOUTME: diffs the children against the package's, checks each child and times both.
# Usage: julia compare.jl <set>...   (default: every set in data/pkg_*.jls)
include(joinpath(@__DIR__, "load.jl"))

const ATOMS_PP = Dict(:A => [:C => 1, :X => 1], :B => [:N => 1], :P => [:C => 1],
                      :Q => [:N => 1, :X => 1])
const ATOMS_BIBI = Dict(:A => [:C => 1], :B => [:N => 1], :P => [:C => 1], :Q => [:N => 1])
const ATOMS_UNI = Dict(:S => [:C => 1], :P => [:C => 1])
const SETS = Dict(
    "init_bibi_pp" => (core = "init_bibi_pp", declared = [:I, :A, :Q], atoms = _ -> ATOMS_PP),
    "dead1" => (core = "", declared = [:I, :A, :Q], atoms = _ -> ATOMS_PP),
    "expand_bibi" => (core = "expand_bibi", declared = [:I, :A, :Q], atoms = _ -> ATOMS_PP),
    "allo_bibi" => (core = "allo_bibi", declared = [:I, :A, :Q], atoms = _ -> ATOMS_BIBI),
    "p1_enum_inh" => (core = "p1_enum_inh", declared = [:I],
                      atoms = m -> :S in m.substrates ? ATOMS_UNI : ATOMS_BIBI))

canon(p) = Slice.Converted(p, 0).mech

"""`f()` and its wall time: the minimum of three runs when one run is under 50 ms (the
package run's protocol)."""
function _testhelper_timed(f)
    t = @elapsed r = f()
    t < 0.05 && for _ in 1:2
        t = min(t, @elapsed f())
    end
    r, t
end

"""Multiset difference `a - b` of two child lists."""
function _testhelper_minus(a, b)
    left = Dict{Slice.Mechanism, Int}()
    foreach(x -> left[x] = get(left, x, 0) + 1, b)
    out = Slice.Mechanism[]
    for x in a
        get(left, x, 0) > 0 ? (left[x] -= 1) : push!(out, x)
    end
    out
end

"""Problems with one child of the slice: not canonical through `raw`, a failed boundary
check or an atom imbalance."""
function _testhelper_child_problems(c, atoms)
    out = String[]
    Slice.canonicalize(Slice.raw(c)[1]...; Slice.raw(c)[2]...)[1] == c ||
        push!(out, "not canonical")
    try Slice.check_boundary(c) catch e; push!(out, "boundary: " * sprint(showerror, e)) end
    try Slice.check_atoms(c, atoms) catch e; push!(out, "atoms: " * sprint(showerror, e)) end
    out
end

const MOVES = (:flip, :split, :dead)

function _testhelper_run(set)
    cfg = SETS[set]
    d = deserialize(joinpath(S, "slice/enum/data/pkg_$set.jls"))
    core = isempty(cfg.core) ? nothing :
           deserialize(joinpath(S, "slice/core/data/$(cfg.core).jls"))
    rows = NamedTuple[]
    nalign = 0
    for (i, r) in enumerate(d.records)
        m = canon(r.parent)
        core === nothing || (nalign += core[d.indices[i]].mech == m)
        fs = (flip = () -> Slice.flip_children(m), split = () -> Slice.split_children(m),
              dead = () -> Slice.dead_end_children(m, cfg.declared))
        row = Dict{Symbol, Any}(:index => d.indices[i], :pkg_count => r.count,
                                :count => try length(Slice.solve(m).indep) catch; -1 end,
                                :allosteric => Slice.is_allosteric(m), :mech => m)
        for mv in MOVES
            pk = getfield(r, mv)
            pkg = canon.(pk.children)
            mine, t, err = try
                (_testhelper_timed(fs[mv])..., "")
            catch e
                (Slice.Mechanism[], NaN, first(sprint(showerror, e), 300))
            end
            probs = [p for c in mine for p in _testhelper_child_problems(c, cfg.atoms(m))]
            counts = mv === :split && isempty(err) ?
                     [try length(Slice.solve(c).indep) catch; -1 end for c in mine] : Int[]
            row[mv] = (pkg = pkg, mine = mine, missing = _testhelper_minus(pkg, mine),
                       extra = _testhelper_minus(mine, pkg), t_pkg = pk.time, t_mine = t,
                       err_pkg = pk.error, err_mine = err, problems = probs,
                       pkg_counts = pk.counts, counts)
        end
        push!(rows, (; (k => v for (k, v) in row)...))
    end
    @printf("\n== %s: %d parents (%d plain, %d allosteric)", set, length(rows),
            count(r -> !r.allosteric, rows), count(r -> r.allosteric, rows))
    core === nothing ? println() : @printf("; %d of %d equal the core's conversion\n",
                                           nalign, length(rows))
    plain = filter(r -> !r.allosteric, rows)
    @printf("  parent counts: plain %d/%d equal the package's; allosteric %d/%d\n",
            count(r -> r.count == r.pkg_count, plain), length(plain),
            count(r -> r.allosteric && r.count == r.pkg_count, rows),
            count(r -> r.allosteric, rows))
    for mv in MOVES
        xs = [r[mv] for r in rows]
        same = count(x -> isempty(x.missing) && isempty(x.extra) && isempty(x.err_pkg) &&
                          isempty(x.err_mine), xs)
        tp = [x.t_pkg for x in xs if isfinite(x.t_pkg)]
        tm = [x.t_mine for x in xs if isfinite(x.t_mine)]
        @printf("  %-5s identical child multisets %d/%d; children pkg %d, mine %d; \
                missing %d, extra %d; errors pkg %d, mine %d; child problems %d\n",
                mv, same, length(xs), sum(length(x.pkg) for x in xs),
                sum(length(x.mine) for x in xs), sum(length(x.missing) for x in xs),
                sum(length(x.extra) for x in xs), count(x -> !isempty(x.err_pkg), xs),
                count(x -> !isempty(x.err_mine), xs), sum(length(x.problems) for x in xs))
        @printf("        time per parent: pkg median %.1f µs (p90 %.1f, total %.3f s); \
                mine median %.1f µs (p90 %.1f, total %.3f s); median ratio %.1f\n",
                1e6 * median(tp), 1e6 * quantile(tp, 0.9), sum(tp), 1e6 * median(tm),
                1e6 * quantile(tm, 0.9), sum(tm),
                median([x.t_pkg / x.t_mine for x in xs if isfinite(x.t_pkg) &&
                        isfinite(x.t_mine) && x.t_mine > 0]))
        if mv === :split
            pairs = [(c, Dict(zip(x.pkg, x.pkg_counts))) for x in xs for c in x.mine]
            known = [(c, d) for (c, d) in pairs if haskey(d, c)]
            mine_counts = Dict(c => n for x in xs for (c, n) in zip(x.mine, x.counts))
            @printf("        split children found in the package's: %d; their counts equal \
                    the package's: plain %d/%d, allosteric %d/%d\n", length(known),
                    count(((c, d),) -> !Slice.is_allosteric(c) && d[c] == mine_counts[c], known),
                    count(((c, d),) -> !Slice.is_allosteric(c), known),
                    count(((c, d),) -> Slice.is_allosteric(c) && d[c] == mine_counts[c], known),
                    count(((c, d),) -> Slice.is_allosteric(c), known))
        end
    end
    serialize(joinpath(S, "slice/enum/data/cmp_$set.jls"), rows)
    rows
end

println("thermo.jl ", _testhelper_src_sha("thermo.jl"), ", enum.jl ",
        _testhelper_src_sha("enum.jl"), ", types.jl ", _testhelper_src_sha("types.jl"))
for set in (isempty(ARGS) ? ["init_bibi_pp", "dead1", "expand_bibi", "allo_bibi",
                             "p1_enum_inh"] : ARGS)
    _testhelper_run(set)
end
