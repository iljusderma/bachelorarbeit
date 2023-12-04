include("main.jl")
using PlotThemes, CSV, Tables
theme(:default)

function animate(snaplen, t0, STATES, α, β, p1)
    ind = 1:snaplen
    SNAPS = zeros((L, Int(t0/snaplen)))
    # create snapshots of the density profile
    for i in range(1, Int(t0/snaplen))
        SNAPS[:, i] = vec(mean(STATES[:, i .* ind], dims=2))
    end
    # animate plots of the densityprofile at each timestep
    anim = @animate for (i, snap) in enumerate(eachcol(SNAPS))
        scatter(snap, ylim=(0, 1), 
        label="α=$α, β=$β, p=$p1", 
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i",
        title="$(i*snaplen) iterations", 
        legend=:topright)
    end
    g = gif(anim, fps=5)
    display(g)
end

# initialize lattice parameters
# lattice size L, injection rate α, ejection rate β, hop rate p
t0 = 400_000 # one time unit includes L updates of the lattice
L = 4000
α = 0.8
β = 0.8
p1 = 1
p2 = 0.3

# perform update
@time begin
#ProfileView.@profview
STATES, CURRENT = simulate(α, β, L, t0, p1, p2)      # HD, LD phase
end

# plot data
cut_STATES = STATES[:, 5_000:end]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))

p = scatter(totaldensity, msw=0, ms=2,  
    label="α=$α, β=$β, p1=$p1, p2=$p2", 
    ylims=[0, 1],
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright)

# P2 = range(p2, 1, 4)[2:end]
# for p2 in P2
#     _STATES, _CURRENT = simulate(α, β, L, t0, p1, p2)
#     _cut_STATES = _STATES[:, 2500:end]
#     _densityprofile = vec(mean(_cut_STATES, dims=2))
#     scatter!(_densityprofile, msw=0, 
#     label="α=$α, β=$β, p1=$p1, p2=$p2")
#     println(p2)
# end

#=
left, right = 1:Int(L/2), Int(L/2)+1:L
plot!(left, fill(α/(p2), Int(L/2)), label="Expected density in left lattice LD")
plot!(right, fill(α, Int(L/2)), label="Expected density in right lattice LD")

plot!(left, fill(1/(p2+1), Int(L/2)), label="Expected density in left lattice MC", ls=:dash)
plot!(right, fill(p2/(p2+1), Int(L/2)), label="Expected density in right lattice MC", ls=:dash)
#h = heatmap(STATES', ylabel="Time t", xlabel="Lattice site i", legend=false, c=:grays, fmt=:png)=#

#animate(500, t0, STATES, α, β, p1)

display("image/png", p) # export as png