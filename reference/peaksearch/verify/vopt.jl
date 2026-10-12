# ABOUTME: Derivative-free optimizers for the verifier searches: differential evolution with
# ABOUTME: restarts on stagnation, and Nelder–Mead.

"""DE (rand/1 and rand-to-best/1, binomial, dithered F, CR in {0.1, 0.9}) with restarts on stagnation,
for tmax seconds."""
function de!(f, lo, hi, rng; NP = 120, tmax = 60.0, seeds = Vector{Vector{Float64}}())
    D = length(lo)
    newpop() = [lo .+ (hi .- lo) .* rand(rng, D) for _ in 1:NP]
    pop = newpop()
    for (i, s) in enumerate(seeds)
        i > NP && break
        pop[i] = clamp.(copy(s), lo, hi)
    end
    fit = [f(x) for x in pop]
    t0 = time()
    trial = zeros(D)
    gens = 0
    restarts = 0
    lastbest = minimum(fit)
    stall = 0
    while time() - t0 < tmax
        gens += 1
        ib = argmin(fit)
        for i in 1:NP
            r1 = rand(rng, 1:NP)
            while r1 == i
                r1 = rand(rng, 1:NP)
            end
            r2 = rand(rng, 1:NP)
            while r2 == i || r2 == r1
                r2 = rand(rng, 1:NP)
            end
            r3 = rand(rng, 1:NP)
            while r3 == i || r3 == r1 || r3 == r2
                r3 = rand(rng, 1:NP)
            end
            F = 0.4 + 0.6 * rand(rng)
            CR = rand(rng) < 0.5 ? 0.1 : 0.9
            tobest = rand(rng) < 0.3
            jr = rand(rng, 1:D)
            for j in 1:D
                if j == jr || rand(rng) < CR
                    base = tobest ? pop[r1][j] + F * (pop[ib][j] - pop[r1][j]) : pop[r1][j]
                    y = base + F * (pop[r2][j] - pop[r3][j])
                    y < lo[j] && (y = lo[j] + rand(rng) * (pop[i][j] - lo[j]))
                    y > hi[j] && (y = hi[j] - rand(rng) * (hi[j] - pop[i][j]))
                    trial[j] = y
                else
                    trial[j] = pop[i][j]
                end
            end
            ft = f(trial)
            if ft <= fit[i]
                pop[i] .= trial
                fit[i] = ft
            end
        end
        b = minimum(fit)
        if b < lastbest - 1e-15
            lastbest = b
            stall = 0
        else
            stall += 1
        end
        spread = maximum(fit) - minimum(fit)
        if stall > 300 || spread < 1e-15
            restarts += 1
            pop = newpop()
            fit = [f(x) for x in pop]
            lastbest = minimum(fit)
            stall = 0
        end
    end
    gens, restarts
end

function nelder_mead(f, x0, step; maxevals = 20_000, tol = 1e-15)
    n = length(x0)
    S = [copy(x0)]
    for i in 1:n
        x = copy(x0)
        x[i] += step
        push!(S, x)
    end
    fs = [f(x) for x in S]
    evals = n + 1
    while evals < maxevals
        o = sortperm(fs)
        S = S[o]
        fs = fs[o]
        abs(fs[end] - fs[1]) < tol && break
        c = sum(S[1:n]) / n
        xr = c + (c - S[end])
        fr = f(xr)
        evals += 1
        if fr < fs[1]
            xe = c + 2 * (c - S[end])
            fe = f(xe)
            evals += 1
            if fe < fr
                S[end], fs[end] = xe, fe
            else
                S[end], fs[end] = xr, fr
            end
        elseif fr < fs[n]
            S[end], fs[end] = xr, fr
        else
            xc = fr < fs[end] ? c + 0.5 * (xr - c) : c + 0.5 * (S[end] - c)
            fc = f(xc)
            evals += 1
            if fc < min(fr, fs[end])
                S[end], fs[end] = xc, fc
            else
                for i in 2:n+1
                    S[i] = S[1] + 0.5 * (S[i] - S[1])
                    fs[i] = f(S[i])
                    evals += 1
                end
            end
        end
    end
    i = argmin(fs)
    S[i], fs[i]
end

