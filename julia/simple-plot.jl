using Plots, Statistics, Random, LaTeXStrings
using PlotThemes
include("main.jl")
theme(:lime)

function plot_densityprofile(densityprofile, L, t, total_density)
    gr()
    p1 = scatter(1:L, densityprofile, 
                    label="$rates", 
                    title="Mean density per lattice site", ms=1, 
                    msw=0, 
                    ylabel=L"\langle \rho_i \rangle", 
                    xlabel="Lattice site "*L"i")
    p2 = scatter(1:(t - L+1), total_density, 
                    label="$rates", title="Time evolution of the density", 
                    ms=1, msw=0)
    plot(p1, p2, layout=(2,1))
end

function state_map(all_states)
    gr()
    heatmap(all_states', 
        xlabel="Lattice site i", ylabel="Time t", 
        cbar=false, title="Movement of the shock over time")
end

function dynamic_plot(all_states)
    @gif for i in 1:t
        plot(all_states[:, i]', title="$i", label=false)
    end every 1
end

L = 300
t = 2*1000
occupied_ratio = 0.5
rates = [0.3, 0.3, 0.75, 0]
flux_steps = 500
state = initialState(L, occupied_ratio)

all_states, FLUX = RANDUpdate(t, state, rates, flux_steps)
# cut data to exclude beginning phase (levelling) which is approx L
println(size(all_states))
cut_states = all_states[:, L:t]
densityprofile = vec(mean(cut_states, dims=2))
total_density = vec(mean(cut_states, dims=1))

# plot data
# plot_densityprofile(densityprofile, L, t, total_density)
state_map(all_states)