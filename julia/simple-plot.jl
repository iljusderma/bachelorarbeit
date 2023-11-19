include("main.jl")
theme(:lime)

# initialize lattice parameters
# lattice size L, injection rate α, ejection rate β, hop rate p
t0 = 100_000 # one time unit includes L updates of the lattice
L = 500
α = 0.8
β = 0.8
p1 = 1
p2 = 1

# perform update
@time begin
#ProfileView.@profview
STATES, curr = simulate(α, β, L, t0, p1, p2)      # HD, LD phase
end

# plot data
cut_STATES = STATES[:, 2500:end]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))
println(curr)


scatter(densityprofile, ms=1, msw=0, 
    label="α=$α, β=$β, p1=$p1, p2=$p2", 
    ylims=[0, 1], 
    title="Totaldensity",
    ylabel=L"\langle \rho_{tot} \rangle", 
    xlabel="Time t", 
    legend=:bottomright)
#hline!([0.7], label="Expected density in HD phase")
#hline!([0.6], label="Expected density in (HD, LD)", ls=:dash)

#heatmap(STATES', ylabel="Time t", xlabel="Lattice site i")