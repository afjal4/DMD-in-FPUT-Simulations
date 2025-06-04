using Plots, LinearAlgebra, RandomizedLinAlg
import DifferentialEquations: ODESolution
include("ode_system.jl")
include("visualisers.jl")

function Ã(svdecomp:: SVD, X_prime:: Matrix{Float64})
    return  svdecomp.U' * X_prime * svdecomp.V * diagm(1 ./ svdecomp.S) #Ut * X' * V * Σ^-1
end

function truncate_SVD(svdecomp:: SVD, r::Int)
    # Truncate the SVD to the first r modes
    U_trunc = svdecomp.U[:, 1:r]
    S_trunc = svdecomp.S[1:r]
    V_trunc = svdecomp.V[:, 1:r]
    return LinearAlgebra.SVD(U_trunc, S_trunc, V_trunc')
end

function mode_decomposition(sol, phases = 1, shift = 0.1, r = 5, use_randsvd = false)
    # Perform Dynamic Mode Decomposition (DMD)
    #=
    args:
    sol: ODESolution object containing the solution to the ODE problem
    phases: Number of phases to decompose into
    shift: % of size of X, to shift to X'
    r: number of modes to extract
    use_randsvd: If true, use randomized SVD, otherwise use standard SVD
    =#

    N = length(sol.u[1]) ÷ 2
    nt =  length(sol.t)
    eigenvalues = []
    modes = []

    phase_size = nt ÷ phases
    Δt = Int(ceil(phase_size * shift))
    X_length = phase_size - Δt

    for i in 1:phases
        # Extract the data for the current phase
        X_phase = hcat([sol.u[t] for t in 
                ((i-1) * phase_size + 1) : (i * phase_size)]...)

        # Split the data into X and X' (overlapping)
        X = X_phase[:, 1:X_length]
        X_prime = X_phase[:, Δt+1:end]

        # Perform SVD on X to find Ã, and it's spectrum
        X_svd = (use_randsvd) ? rsvd(X, r) : truncate_SVD(svd(X), r)
        Λ, W = eigen(Ã(X_svd, X_prime))
        modes_from_phase = X_svd.U * W  # Modes are U * W
        
        push!(eigenvalues, Λ)
        push!(modes, modes_from_phase)
    end
    return eigenvalues, modes
end

N = 500
t = 10000
α = 0.5
β = 0.75

phases = 100
shift = 0.1 # size /1 of shift from X to X'

sol = solve(ODE_problem(N, α, β, t), 
            KenCarp47(linsolve = KrylovJL_GMRES()),
            saveat=1, 
            reltol=1e-6, abstol=1e-6)

evals_progression, modes_progression = mode_decomposition(sol, phases, shift)

eigenvalues_plot(evals_progression, show = true)