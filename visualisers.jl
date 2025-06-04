using Plots
import DifferentialEquations: ODESolution

function displacement_gif(sol::ODESolution)
    N = length(sol.u[1]) ÷ 2
    anim = @animate for i in 1:length(sol.t)
        q = sol.u[i][1:N]
        # Create array including fixed endpoints
        q_full = [0.0; q; 0.0]  # [q₀, q₁, ..., qₙ, qₙ₊₁]
        plot(0:N+1, q_full, ylim=(-0.5, 0.5), lw=3, marker=:circle, ms=6,
             xlabel="Position in chain", ylabel="Displacement",
             title="FPUT Chain (N=$(N)) with Fixed Endpoints at t=$(round(sol.t[i], digits=2))",
             legend=false)
        # Mark the fixed endpoints
        scatter!([0, N+1], [0, 0], color=:red, marker=:circle, ms=8, label="Fixed Endpoints")
        hline!([0], lw=1, linestyle=:dash, color=:black)
    end
    return anim
end

function displacement_heatmap(sol::ODESolution; show::Bool=false)
    N = length(sol.u[1]) ÷ 2
    # Create a matrix of displacements over time
    # Each row is a time point, each column is a particle position
    displacement_matrix = zeros(length(sol.t), N)
    for (i, t) in enumerate(sol.t)
        displacement_matrix[i, :] = sol.u[i][1:N]
    end
    
    hmap = heatmap(displacement_matrix, 
                  xlabel="Position in chain", 
                  ylabel="Time step",
                  title="Displacement Heatmap",
                  aspect_ratio=:auto,
                  colorbar_title="Displacement")
    if show
        display(plt)
    end
    return hmap
end

function hamiltonian(sol::ODESolution, α, β)
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

function hamiltonian_plot(sols::Vector{ODESolution}, α, β; show::Bool=false)
    plt = plot()
    plt.title = "Total Energy of FPUT Chain Over Time"

    for (i, sol) in enumerate(sols)
        plot!(plt, sol.t, hamiltonian(sol, α, β), lw=2, label="Solution $i")
    end
    xlabel!(plt, "Time")
    ylabel!(plt, "Total Energy")
    title!(plt, "Total Energy of FPUT Chain Over Time")
    if show
        display(plt)
    end
    return plt
end

function eigenvalues_plot(eigenvalues::Vector; show::Bool=false)
    plt = plot()
    phases = length(eigenvalues)
    
    
    colours = cgrad(:bluesreds, phases)

    for (i, Λ) in enumerate(eigenvalues)
        plot!(plt, real.(Λ), imag.(Λ), 
        seriestype=:scatter, 
        color= colours[i], 
        markersize=2,
        legend = false)
    end
    xlabel!(plt, "Real Part")
    ylabel!(plt, "Imaginary Part")
    title!(plt, "Eigenvalues of FPUT Chain Over Time")
    if show
        display(plt)
    end
    return plt
end

function transformed_eigenvalues_plot(eigenvalues::Vector, Δt; show::Bool=false)
    plt = plot()
    phases = length(eigenvalues)
    
    colours = cgrad(:bluesreds, phases)

    transformed_evals = [ln.(evals) ./ Δt for evals in eigenvalues]
    for (i, Λ) in enumerate(transformed_evals)
        plot!(plt, real.(Λ), imag.(Λ), 
        seriestype=:scatter, 
        color= colours[i], 
        markersize=2,
        legend = false)
    end
    xlabel!(plt, "Real Part; Growth Rate")
    ylabel!(plt, "Imaginary Part; Frequency")
    title!(plt, "Transformed Eigenvalues of FPUT Chain Over Time")
    if show
        display(plt)
    end
    return plt
end

    