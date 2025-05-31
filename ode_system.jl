# FPUT RHS function: du/dt = f(u,t)
function fput_rhs!(du, u, p, t)
    N, α, β = p               # number of particles, nonlinearity coefficients
    q = u[1:N]             # positions
    pvec = u[N+1:2N]       # momenta

    dq = pvec              # dq/dt = p
    dp = zeros(N)

    # Fixed boundary conditions: q_0 = q_{N+1} = 0
    q_left = 0.0
    q_right = 0.0

    for j in 1:N
        # Use neighbors or boundaries
        q_prev = (j == 1) ? q_left : q[j-1]
        q_next = (j == N) ? q_right : q[j+1]

        dp[j] = (q_next - 2q[j] + q_prev) 
                + α * ((q_next - q[j])^2 - (q[j] - q_prev)^2)
                + β * ((q_next - q[j])^3 - (q[j] - q_prev)^3)
    end

    du[1:N] .= dq
    du[N+1:2N] .= dp
end

#Wrapper function to set up the ODE problem
function ODE_problem(N, α, β, t)
    p = (N, α, β)
    q0 = zeros(N)
    p0 = zeros(N)

    #Initial Preturbation
    q0 = [0.1 * sin(2*pi*i/N) for i in 1:N]

    u0 = [q0; p0]
    tspan = (0.0, t)
    prob = ODEProblem(fput_rhs!, u0, tspan, p)
    return prob
end
