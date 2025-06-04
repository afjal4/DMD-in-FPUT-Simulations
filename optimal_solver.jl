#Optimised Solver using Symplectic  for FPUT simulations
#Using RandSVD to find dominant modes

using DifferentialEquations, Plots
include("visualisers.jl")
include("ode_system.jl")

function A(v)
    return circshift(v, -1) - 2 * v + circshift(v, 1)
end

function B(v)
    return circshift(v, -1) - circshift(v, 1)
end

function jacobian_vector_prod(v, α)
    N = length(u) ÷ 2
    Jv = similar(u)

    # dF/dq = 1
    Jv[1:N] .= v[N+1:2N]

    # dF/dp = f'(q) = A(1 + Bαq) + BAαq, f(q) = (Aq).(1 + αBq) 
    q = v[1:N]
    Jv[N+1:2N] .= A(1 .+ B(α * q)) + α * B(A(q))
    return Jv
end

function jvp!(Jv, v, u, p, t)
    # Compute the Jacobian-vector product for the FPUT system
    Jv .= jacobian_vector_prod(u, p[2])
end