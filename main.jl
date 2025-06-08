#Naive Simulation for small N
#SVD to find modes ##WIP

using DifferentialEquations, Plots
include("visualisers.jl")
include("ode_system.jl")

# Parameters
N = 100             # number of particles
α = 1              # nonlinearity coefficients:
β = 0.5             # 
t = 1000.0            # time span



prob = ODE_problem(N, α, β, t)
sol = solve(prob, KenCarp47(linsolve = KrylovJL_GMRES()), dt=0.01)
#sol2 = solve(prob, KenCarp47(linsolve = KrylovJL_GMRES()), reltol=1e-6, abstol=1e-6)

#GIF
#=
anim = displacement_gif(sol)
display(gif(anim, fps=30))
=#

#Plots
#display(displacement_heatmap(sol1))
display(hamiltonian_plot(sol, α, β))
#display(total_energy_plot([sol1,sol2], α, β))
