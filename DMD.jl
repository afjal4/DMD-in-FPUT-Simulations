using Plots, RandomizedLinAlg
import DifferentialEquations: ODESolution

function randsvd(X, X', r)
    # Perform randomized SVD
    #=
    args:
    X: Matrix of data
    X': Matrix of shifted data
    r: Target rank for the SVD
    =#

    #Truncates the matrices to the minimum size if they are not equally sized
    size = min(size(X,2), size(X',2))
    X = @view X[:, 1:size]
    X'= @view X'[:, 1:size]
    
end




function mode_decomposition(sol, phases, shift, r = -1)
    # Perform Dynamic Mode Decomposition (DMD)
    #=
    args:
    sol: ODESolution object containing the solution to the ODE problem
    phases: Number of phases to decompose into
    shift: % of size of X, to shift to X'
    =#

    # Method of SVD in DMD determined by whether r is specified 
    use_randsvd = (r<=0) ? false : true

    N = length(sol.u[1]) ÷ 2
    nt = length(sol.t)
    results = []

    phase_size = nt ÷ phases
    Δt = floor(phase_size * shift)
    X_length = phase_size - Δt

    for i in 1:phases
        X_phase = hcat([sol.u[i][1:N] for j in 
                (i * phase_size) : ((i+1) * phase_size - 1)])

        X = X_phase[:, 1:X_length]
        X' = X_phase[:, Δt+1:end]

        result = (use_randsvd) ? 
            randsvd(X, X', r) : svd(X, X')

        push!(results, result)
    end
    return results
end

N = 1000
t = 1000
α = 0.5
β = 0.5

sol = solve(ODE_problem(N, α, β, t), 
            KenCarp47(linsolve = KrylovJL_GMRES()),    
            saveat=1, 
            reltol=1e-6, abstol=1e-6)

phases = 5
shift = 0.1 # size /1 of shift from X to X'


X = hcat()