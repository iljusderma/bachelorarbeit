include("main.jl")
using ProfileView
theme(:lime)

# initialize lattice parameters
# lattice size L, injection rate α, ejection rate β, hop rate p
t0 = 100_000 # one time unit includes L updates of the lattice
L = 500
α = 0.3
β = 0.8
p1 = 1
p2 = 1

# perform update
@time begin
#ProfileView.@profview
STATES, curr = simulate(α, β, L, t0, p1, p2)
end
# plot data
cut_STATES = STATES[:, 2500:t0]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))


scatter(totaldensity, ms=1, 
    label="α=$α, β=$β, p=$p", 
    ylims=[0, 1], 
    title="Density profile",
    ylabel=L"\langle \rho \rangle", 
    xlabel="Lattice site i")
hline!([0.7], label="Expected density in HD phase")
hline!([0.3], label="Expected density in LD phase")

#heatmap(STATES', ylabel="Time t", xlabel="Lattice site i")