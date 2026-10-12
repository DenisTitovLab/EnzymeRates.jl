# ABOUTME: Prints worked examples of each census class from saved runs: key duplicates with the
# ABOUTME: candidate they repeat, children with their parent, W, P by sub-class, and Xlaw names.

include(joinpath(@__DIR__, "census_lib.jl"))
include(joinpath(@__DIR__, "summarize.jl"))

const NEX = parse(Int, ARGS[1])
for name in Symbol.(ARGS[2:end])
    data = deserialize(joinpath(@__DIR__, "..", "data", "$(name).jls"))
    rounds = data.rounds
    label!(rounds)
    byid = Dict(r[:id] => r for recs in rounds for r in recs)
    first_with_key = Dict{String, Any}()
    for recs in rounds, r in recs
        r[:fitted] && haskey(r, :key) && get!(first_with_key, r[:key], r)
    end
    fit = [r for recs in rounds for r in recs if r[:fitted]]
    println("="^100, "\n", name)
    show(r) = @sprintf("id %s, %d package params (rank %d), %d slice params (rank %d), D terms %d\n        %s",
                       r[:id], r[:n_params], get(r, :pkg_rank_used, -1), get(r, :n_slice, -1),
                       get(r, :rank, -1), get(r, :nD, -1), get(r, :text, ""))
    # D: who they repeat, by kind of pair
    D = [r for r in fit if :D in r[:labels]]
    kinds = Dict{String, Int}()
    for r in D
        o = first_with_key[r[:key]]
        k = string(o[:id][1] == r[:id][1] ? "same round" : "earlier round",
                   o[:id] == r[:parent] ? ", the parent" :
                   (r[:parent] != 0 && get(o, :parent, 0) == r[:parent]) ? ", a sibling" : "",
                   get(r, :allo, false) ? ", allosteric" : ", plain",
                   get(o, :allo, false) ? " = allosteric" : " = plain")
        kinds[k] = get(kinds, k, 0) + 1
    end
    println("D (key duplicates of an earlier fitted candidate): ", length(D), " by pair kind: ", kinds)
    for r in D[1:min(NEX, end)]
        println("  dup  ", show(r), "\n  of   ", show(first_with_key[r[:key]]))
    end
    Z = [r for r in fit if :Z in r[:labels]]
    println("Z (same key as the parent): ", length(Z))
    for r in Z[1:min(NEX, end)]
        println("  child  ", show(r), "\n  parent ", show(byid[r[:parent]]))
    end
    for l in (:W_admit, :W_split, :W_trap, :W_bottomless, :W_D0_face, :W_v0, :W_re_turnover,
              :W_boundary, :W_ctor, :W_law_err)
        rs = [r for r in fit if l in r[:labels]]
        isempty(rs) && continue
        println(l, ": ", length(rs))
        for r in rs[1:min(NEX, end)]
            msg = l === :W_admit ? r[:admit] : l === :W_boundary ? r[:boundary] :
                  l === :W_D0_face ? string(r[:dead]) : l === :W_v0 ? "N ≡ 0" :
                  l === :W_re_turnover ? string(r[:re_turn]) :
                  string(get(r, :faces, ""), get(r, :conv_err, ""), get(r, :law_err, ""))
            println("  ", show(r), "\n        verdict: ", first(msg, 400))
        end
    end
    P = [r for r in fit if :P in r[:labels] || :Ptoday in r[:labels]]
    subs = Dict{String, Vector{Any}}()
    for r in P
        k = join(sort!(string.(collect(keys(get(r, :explain, Dict{Symbol, Int}()))))), "+")
        isempty(k) && (k = "package only (slice full rank)")
        push!(get!(subs, k, Any[]), r)
    end
    println("P or Ptoday by slice attribution: ", Dict(k => length(v) for (k, v) in subs))
    for (k, rs) in sort!(collect(subs); by = first)
        println("  [$k]")
        for r in rs[1:min(NEX, end)]
            println("  ", show(r), "  ", get(r, :explain, ""))
        end
    end
    X = [r for r in fit if :Xlaw in r[:labels]]
    println("Xlaw (package law rank differs from slice law rank): ", length(X),
            "; package count - slice count: ",
            sort!(collect(Dict(d => count(r -> r[:n_params] - r[:n_slice] == d, X)
                               for d in unique(r[:n_params] - r[:n_slice] for r in X)))))
    for r in X[1:min(NEX, end)]
        em = ER.compile_mechanism(r[:mech])
        _, indep = ER._dependent_param_exprs(typeof(em))
        sl = slice_fitted(r[:smech], S.solve(r[:smech]))
        norm(x) = replace(string(x), "_res_+" => "_res_plus", "_-" => "_minus")
        pk, sk = Set(norm.(indep)), Set(string.(sl))
        println("  ", show(r))
        println("        package only: ", sort!(collect(setdiff(pk, sk))),
                "; slice only: ", sort!(collect(setdiff(sk, pk))))
    end
    N = [r for r in fit if :N in r[:labels]]
    println("N: ", length(N))
    for r in N[1:min(NEX, end)]
        println("  ", show(r), " fwd ", r[:nosat_fwd], " rev ", r[:nosat_rev])
    end
    flush(stdout)
end
