using Plots
import DifferentialEquations: ODESolution

function total_energy(sol::ODESolution, α, β)
    N = length(sols.u) ÷ 2
    nt = zeros(length(sols.t))

    energy = zeros(nt)

    # Calculate total energy at each time step
    for j in 1:nt
        u = sol.u[j]
        q = @view u[1:N]
        p = @view u[N+1:end]
        kinetic_energy = sum(p.^2) / 2
        #Summation over differences between each particle
        potential = 0.0
        for i in 1:(N-1)
            dq = q[i] - q[i+1]
            potential += 0.5 * dq^2 + (α/3) * dq^3 + (β/4) * dq^4
        end
        energy[j] = kinetic_energy + potential_energy
    end
    return energy
end

function total_energy_plot(sols::Vector{ODESolution}, α, β)
    plt = plot()
    plt.title = "Total Energy of FPUT Chain Over Time"

    for (i, sol) in enumerate(sols)
        plot!(plt, sol.t, total_energy(sol, α, β), lw=2, label="Solution $i")
    end
    xlabel!(plt, "Time")
    ylabel!(plt, "Total Energy")
    title!(plt, "Total Energy of FPUT Chain Over Time")
    display(plt)
end
