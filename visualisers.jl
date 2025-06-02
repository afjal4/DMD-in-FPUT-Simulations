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

function displacement_heatmap(sol::ODESolution)
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

    return hmap
end

