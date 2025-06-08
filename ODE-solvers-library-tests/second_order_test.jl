using DifferentialEquations
using Base.Threads

# Parameters
N = 1_000_000                # Number of particles
α = 0.25                     # Nonlinearity parameter

# Threaded acceleration function
function fput_accel_threaded!(d2x, x, dx, p, t)
    α = p
    N = length(x)

    @threads for j in 2:N-1
        lin = x[j-1] + x[j+1] - 2x[j]
        nonlin = (x[j+1] - x[j])^3 - (x[j] - x[j-1])^3
        d2x[j] = lin + α * nonlin
    end

    # Boundary conditions (fixed ends)
    d2x[1] = 0.0
    d2x[N] = 0.0
end

# Initial conditions
x0 = [sin(2 * π * i / N) for i in 1:N]  # Initial displacement
v0 = zeros(N)

# Time span
tspan = (0.0, 10000.0)

# Define the problem
prob = SecondOrderODEProblem(fput_accel_threaded!, x0, v0, tspan, α)

# Solve using Velocity Verlet (symplectic integrator)
sol = solve(prob, VelocityVerlet(); dt=0.5, saveat=10.0)


println("Final displacement at center: ", sol[end][N ÷ 2])
