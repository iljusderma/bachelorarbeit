include("main.jl")
theme(:lime)

# initialize lattice parameters
# lattice size L, injection rate α, ejection rate β, hop rate p
t0 = 100_000 # one time unit includes L updates of the lattice
L = 500
α = 0.3
β = 0.6
p1 = 1
p2 = 1

# perform update
@time begin
#ProfileView.@profview
STATES1, curr1 = simulate(α, β, L, t0, p1, 0.4)      # HD, LD phase
end

# plot data
cut_STATES = STATES1[:, 2500:end]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity1 = vec(mean(STATES1, dims=1))
println(mean(cut_STATES))


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