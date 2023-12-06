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
t0 = 50_000 # one time unit includes L updates of the lattice
L = 200
α = 0.8
β = 0.8
p1 = 1
p2 = 0.3

# perform update
@time begin
STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
end

# plot data
cut_STATES = STATES[:, 3_000:end]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))

p = scatter(densityprofile, msw=0, ms=2,  
    label=L"t=10^4", 
    ylims=[0, 1],
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright)

#animate(500, t0, STATES, α, β, p1)

display("image/png", p) # export as png