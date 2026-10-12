# decode cost per distinct NTuple length
dec1(@nospecialize(code)) = Int[code[i] for i in 1:length(code)]
dec2(@nospecialize(code)) = Int[code...]
dec3(@nospecialize(code)) = collect(Int, code)
for (name, f) in (("index", dec1), ("splat", dec2), ("collect", dec3))
    ts = Float64[]
    for n in (5000, 5001, 5002, 8000, 8001)
        t = Tuple(rand(1:9, n))
        push!(ts, @elapsed f(t))
    end
    t = Tuple(rand(1:9, 5000))
    println(name, ": first calls ", round.(ts; digits = 4), " repeat ", @elapsed f(t))
end
