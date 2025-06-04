using DifferentialEquations, Plots
include("optimal_solver.jl")

# FPUT RHS function: du/dt = f(u,t)
function fput_rhs!(du, u, p, t)
    N, α, β = p               # number of particles, nonlinearity coefficients
    q = u[1:N]             # positions
    pvec = u[N+1:2N]       # momenta

    # du/dt = [dq/dt; dp/dt]
    # 1) dq/dt = p  
    dq = pvec              # dq/dt = p
    du[1:N] .= dq

    # 2) dp/dt = f(q)

    dp = similar(q)  # dp/dt = f(q)

    # Fixed boundary conditions: q_0 = q_{N+1} = 0
    q_left = 0.0
    q_right = 0.0

    @inbounds for j in 1:N
        # Use neighbors or boundaries
        q_prev = (j == 1) ? q_left : q[j-1]
        q_next = (j == N) ? q_right : q[j+1]

        dp[j] = (q_next - 2q[j] + q_prev) 
                + α * ((q_next - q[j])^2 - (q[j] - q_prev)^2)
                + β * ((q_next - q[j])^3 - (q[j] - q_prev)^3)
    end

    du[N+1:2N] .= dp
end

#Wrapper to set up the ODE problem
function ODE_problem(N, α, β, t; use_jvp = false)
    p = (N, α, β)
    q0 = zeros(N)
    p0 = zeros(N)

    #Initial Preturbation
    q0 = [0.1 * sin(4*pi*i/N) for i in 1:N]

    u0 = [q0; p0]
    tspan = (0.0, t)

    f = (use_jvp) ? ODEFunction(fput_rhs!, jvp = jvp!) :
                  ODEFunction(fput_rhs!)

    prob = ODEProblem(f, u0, tspan, p)
    return prob
end

# Second order ODE problem for FPUT
function rhs_2nd_order!(d2q, q, p, t)
    N, α, β = p
    q_left, q_right = 0.0, 0.0

    dp = similar(q)  # d2q/d2t = f(q)

    for j in 1:N
        q_prev = (j == 1) ? q_left : q[j-1]
        q_next = (j == N) ? q_right : q[j+1]
        dp[j] = (q_next - 2q[j] + q_prev) +
                 α * ((q_next - q[j])^2 - (q[j] - q_prev)^2) +
                 β * ((q_next - q[j])^3 - (q[j] - q_prev)^3)
    end

    d2q[1:N] .= dp
    return d2q
end

function second_order_ODE_problem(N, α, β, t)
    p = (N, α, β)
 
    p0 = zeros(N)
    q0 = [0.1 * sin(2*pi*i/N) for i in 1:N] #Initial Preturbation

    tspan = (0.0, t)
    prob = SecondOrderODEProblem(rhs_2nd_order!, q0, p0, tspan, p)
    return prob
end