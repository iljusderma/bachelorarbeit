using Plots, Statistics, Random, LaTeXStrings
using PlotThemes
theme(:juno)
include("main.jl")

L = 1000
t = 100*1000
occupied_ratio = 0.5
rates = [0.2, 0.7, 1, 0]
flux_steps = 500
state = initialState(L, occupied_ratio)

all_states, FLUX = SLUpdate(t, state, rates, flux_steps)
# cut data to exclude beginning phase (levelling) which is approx L
cut_states = all_states[:, L:t]
densityprofile = vec(mean(cut_states, dims=2))
total_density = mean(cut_states, dims=1)

# plot data
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
                
#print
#flux = mean(FLUX)
#exp = rates[1]*(1-rates[1])
#"flux=$flux, exp.=$exp"