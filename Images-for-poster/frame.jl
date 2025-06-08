#making image for poster

# Soliton-like
q = [exp(-(i-N/2)^2/8) * sin(2π*i/N) for i in 1:32]

# Plotting
plt = plot(layout=(2,1), size=(800,600), dpi=300)

plot(q, 
      line=:stem, marker=:circle, color=:royalblue,
      xlabel="Oscillator Index", ylabel="Displacement qᵢ",
      title="FPUT Chain State (N=$N)", 
      label="Displacement")