#Naive Simulation for small N
#SVD to find modes

using DifferentialEquations, Plots
include("visualisers.jl")
include("ode_system.jl")
gr()

# Parameters
N = 1000               # number of particles
α = 0.0              # nonlinearity coefficients:
β = 0.25             # 
t = 1000.0            # time span

prob = ODE_problem(N, α, β, t)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# Display the animation
#animation = displacement_gif(sol)
#display(gif(animation, fps=30))

display(displacement_heatmap(sol))
