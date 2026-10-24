# ABOUTME: Prints the thread and GC-thread counts and the CPU this Julia runs with.
println((VERSION, Threads.nthreads(), Threads.ngcthreads(), Sys.CPU_NAME, Base.JLOptions().opt_level))
