include("main.jl")

# initialize lattice parameters
# lattice size L, injection rate α, ejection rate β, hop rate p
t0 = 100_000 # one time unit includes L updates of the lattice
L = 500
α = 0.3
β = 0.8
p1 = 1
p2 = 0.7

# perform update
@time begin
#ProfileView.@profview
STATES, CURRENT = simulate(α, β, L, t0, p1, p2)      # HD, LD phase
end

# plot data
cut_STATES = STATES[:, 2500:end]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))
println(mean(CURRENT[10_000:end]))

p = scatter(densityprofile, ms=1, msw=0, 
    label="α=$α, β=$β, p1=$p1, p2=$p2", 
    ylims=[0, 1], 
    title="Densityprofile",
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright,
    size=(900, 600))
#=
P2 = range(0.4, p2, 6)[1:end-1]
for p2 in P2
    _STATES, _CURRENT = simulate(α, β, L, t0, p1, p2)
    _cut_STATES = _STATES[:, 2500:end]
    _densityprofile = vec(mean(_cut_STATES, dims=2))
    scatter!(_densityprofile, ms=1, msw=0, 
    label="α=$α, β=$β, p1=$p1, p2=$p2")
    println(p2)
end
display(p)
=#
left, right = 1:Int(L/2), Int(L/2)+1:L
plot!(left, fill(α/(p2), Int(L/2)), label="Expected density in left lattice LD")
plot!(right, fill(α, Int(L/2)), label="Expected density in right lattice LD")

plot!(left, fill(1/(p2+1), Int(L/2)), label="Expected density in left lattice MC", ls=:dash)
plot!(right, fill(p2/(p2+1), Int(L/2)), label="Expected density in right lattice MC", ls=:dash)
#heatmap(STATES', ylabel="Time t", xlabel="Lattice site i")