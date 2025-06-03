#Naive Simulation for small N
#SVD to find modes ##WIP

using DifferentialEquations, Plots
include("visualisers.jl")
include("ode_system.jl")

# Parameters
N = 100             # number of particles
α = 0.5              # nonlinearity coefficients:
β = 0.5             # 
t = 100.0            # time span

prob = ODE_problem(N, α, β, t)
sol1 = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
#sol2 = solve(prob, KenCarp47(linsolve = KrylovJL_GMRES()), reltol=1e-6, abstol=1e-6)

#GIF
#=
anim = displacement_gif(sol)
display(gif(anim, fps=30))
=#

#Plots
display(displacement_heatmap(sol1))
#display(total_energy_plot(sol2, α, β))
#display(total_energy_plot([sol1,sol2], α, β))
