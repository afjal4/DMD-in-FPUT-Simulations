#Optimised Solver using Symplectic  for FPUT simulations
#Using RandSVD to find dominant modes

using DifferentialEquations, Plots
include("visualisers.jl")
include("benchmarking.jl")
include("ode_system.jl")

# Parameters
N = 100             # number of particles
α = 0.5              # nonlinearity coefficients:
β = 0.5             # 
t = 100.0            # time span

solver = KenCarp47(linsolve = KrylovJL_GMRES())

prob = ODE_problem(N, α, β, t)
prob2 = second_order_ODE_problem(N, α, β, t)

#sol = solve(prob, solver, abstol=1e-8, reltol=1e-6)
#sol = solve(prob2, VelocityVerlet(), dt = 0.01)
sol = solve(prob)

println(1)

display(displacement_heatmap(sol))
display(total_energy_plot(sol, α, β))