using Random, Plots, Statistics, PlotThemes, LaTeXStrings, ProfileView
theme(:lime)

function initialize_state()
    # initialize state
    state = zeros(L)
    state[1:Int(floor(L*0.4))] .= 1
    return shuffle(state)
end

function update(state, hop_counter)
    L = length(state)
    site = rand(0:L)
    if site == 0
        if state[1] == 0 && rand() < α
            # injection
            state[1] += 1
        end
    elseif site == L 
        if state[L] == 1 && rand() < β
            # ejection
            state[L] -= 1
        end
    else
        if state[site] == 1 && state[site + 1] == 0 && rand() < p
            # hop to right
            state[site + 1] += 1
            state[site] -= 1
            if site == L/2
                hop_counter += 1
            end
        end
    end
    return state, hop_counter
end

function simulate()
    # initialize state
    state = initialize_state()
    # save state history in all_states every L-th time step
    STATES = zeros(L, t0)
    snapshot = zeros(L, L)
    # insert Current measurement
    hop_counter = 0
    #CURRENT_SNAPS = zeros(t0)
    # perform L*t0 update steps
    for t in 1:(t0*L)
        state, hop_counter = update(state, hop_counter)
        # save snapshot
        if t%L == 0
            STATES[:, div(t, L)] = vec(mean(snapshot, dims=2))
            #CURRENT_SNAPS[div(t, L)] = hop_counter/L
            #hop_counter = 0
        else
            snapshot[:, t%L] = state 
        end
        # count hops

    end
    return STATES, hop_counter/10_000
end
# initialize lattice parameters
# lattice size L, injection rate α, ejection rate β, hop rate p
t0 = 10_000 # one time unit includes L updates of the lattice
L = 200
α = 0.25
β = 0.75
p = 1

# perform update
@time begin
@profview
    STATES, curr = simulate()
end

# plot data
gr()
cut_STATES = STATES[:, 2500:t0]
densityprofile = vec(mean(cut_STATES, dims=2))
totaldensity = vec(mean(STATES, dims=1))
println(curr)
scatter(densityprofile, ms=1, 
    label="α=$α, β=$β, p=$p", 
    ylims=[0, 1], 
    title="Density profile",
    ylabel=L"\langle \rho \rangle", 
    xlabel="Lattice site i")
# hline!([0.8], label="Expected density")
# plot!(zeros(t0) .+ 0.3)