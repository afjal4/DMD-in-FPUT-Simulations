const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences1D

@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3);
else
    @init_parallel_stencil(Threads, Float64, 3);
end

@parallel function test_d!(out, x)
    @inn(out) = @inn(x)
    return
end

x = Data.Array([1,2,3,4,5])
out = @zeros(5)

test_d!(out, x) #!!!!why does this not work?
println(out)