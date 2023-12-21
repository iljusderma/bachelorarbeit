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
t0 = 100_000 # one time unit includes L updates of the lattice
L = 500
α = 0.4
β = 0.8
p1 = 1
p2 = 1

# perform update
@time begin
STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
end

# plot data
cut_STATES = STATES
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))

plot1 = scatter(densityprofile, msw=0, ms=2,  
    label="α=$α, β=$β", title="a)", 
    ylims=[0, 1],
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright)

α=0.8
@time begin
STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
end
# plot data
cut_STATES = STATES
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))

plot2 = scatter(densityprofile, msw=0, ms=2,  
    label="α=$α, β=$β", title="b)", 
    ylims=[0, 1],
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright)

β=0.4
@time begin
STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
end
# plot data
cut_STATES = STATES
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))

plot3 = scatter(densityprofile, msw=0, ms=2,  
    label="α=$α, β=$β", title="c)", 
    ylims=[0, 1],
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright)

α=0.3
β=0.3
STATES, CURRENT = simulate(α, β, L, 5_000, p1, p2)
densityprofile = vec(mean(STATES, dims=2))

plot4 = scatter(densityprofile, msw=0, ms=2,  
    label="α=$α, β=$β", title="d)", 
    ylims=[0, 1],
    ylabel=L"\langle \rho_i \rangle", 
    xlabel="Lattice site i", legend=:topright)


#animate(500, t0, STATES, α, β, p1)

#plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{p_2+1}")
#plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{p_2}{p_2+1}")
#println(mean(CURRENT))
plot(plot1, plot2, plot3, plot4, layout = 4)
savefig("plot.pdf")