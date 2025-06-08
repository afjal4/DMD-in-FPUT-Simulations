using DifferentialEquations, LinearAlgebra, SparseArrays, BenchmarkTools, LinearSolve
include("optimal_solver.jl")

function benchmark_fput_solvers(N::Int; α=0.25, T=10.0, dt=0.01)
    println("Setting up FPUT system with N = $N, α = $α")

    # Initial conditions
    x0 = zeros(N)
    v0 = zeros(N)
    x0[N ÷ 2] = 1.0
    u0 = [x0; v0]  # Full state vector: [x; v]

    # Tridiagonal Jacobian approximation (for preconditioner)
    function J_nonlin_mv!(out, vx)
        out .= 0.0
        for j in 2:N-1
            out[j] = (vx[j-1] + vx[j+1] - 2vx[j]) +
                     α * (3 * (vx[j+1] - vx[j])^2 - 3 * (vx[j] - vx[j-1])^2)
        end
    end

    # Jacobian-vector product (block structured)
    function Jv!(out, v)
        vx = view(v, 1:N)
        vv = view(v, N+1 : 2*N)
        ox = view(out, 1:N)
        ov = view(out, N+1: 2*N)
        copy!(ox, vv)
        J_nonlin_mv!(ov, vx)
    end

    # ODE RHS
    function f!(du, u, p, t)
        x = view(u, 1:N)
        v = view(u, N+1:2*N)
        dx = view(du, 1:N)
        dv = view(du, N+1:2*N)

        copy!(dx, v)
        dv .= 0.0
        for j in 2:N-1
            dv[j] = x[j-1] + x[j+1] - 2x[j] +
                    α * ((x[j+1] - x[j])^3 - (x[j] - x[j-1])^3)
        end
    end

    # ODE Problem
    tspan = (0.0, T)
    prob = ODEProblem(f!, u0, tspan)


    # Solvers to benchmark
    solvers = Dict(
        "KenCarp47(linsolve = KrylovJL_GMRES())" => () -> solve(prob, KenCarp47(linsolve = KrylovJL_GMRES()), dt=dt),
        #"Tsit5"     => () -> solve(prob, Tsit5(), dt=dt),
        #"Rodas5"    => () -> solve(prob, Rodas5(), dt=dt),
        "TRBDF2"    => () -> solve(prob, TRBDF2(), dt=dt)
    )

    # Benchmark
    results = Dict()
    println("Benchmarking...")
    for (name, solverfun) in solvers
        print("  → $name: "); flush(stdout)
        t = @belapsed $solverfun()
        println("$(round(t, digits=3)) s")
        results[name] = t
    end

    return results
end