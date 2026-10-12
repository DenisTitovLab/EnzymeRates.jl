# ABOUTME: Runs the verifier scripts one after another inside a single Julia process (one lock
# ABOUTME: acquisition), each in its own module with its own log; a failing part does not stop the rest.

const DIR = @__DIR__
const PARTS = [(:CheckPoints, "check_points.jl", String[], "check_points.log"),
               (:FirstPassage, "fp.jl", String[], "fp.log"),
               (:SearchT, "vsearch.jl", [get(ENV, "VBUDGET", "1500")], "vsearch.log"),
               (:Extra, "vextra.jl", [get(ENV, "VTB", "60")], "vextra.log")]
const ONLY = isempty(ARGS) ? Symbol[] : Symbol.(split(ARGS[1], ","))

status(msg) = open(io -> println(io, msg), joinpath(DIR, "runall.status"), "a")

for (name, file, args, log) in PARTS
    isempty(ONLY) || name in ONLY || continue
    status("$(name) start $(time())")
    open(joinpath(DIR, log), "w") do io
        redirect_stdout(io) do
            redirect_stderr(io) do
                empty!(ARGS)
                append!(ARGS, args)
                t = time()
                path = joinpath(DIR, file)
                ex = :(module $name
                       include($path)
                       end)
                try
                    Core.eval(Main, ex)
                catch err
                    println("PART FAILED: ", sprint(showerror, err, catch_backtrace()))
                end
                println("part wall time ", round(time() - t; digits = 1), " s")
                flush(stdout)
            end
        end
    end
    status("$(name) done $(time())")
end
status("all done $(time())")
