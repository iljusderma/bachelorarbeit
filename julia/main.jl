using Random, Plots, Statistics, PlotThemes, LaTeXStrings, ProfileView
theme(:lime)

function initialize_state(L)
    # initialize state
    state = zeros(L)
    state[1:Int(floor(L*0.5))] .= 1
    return shuffle(state)
end

function update(state, hop_counter, α, β, p)
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

function simulate(α, β, L, t0, p=1)
    # initialize state
    state = initialize_state(L)
    # save state history in all_states every L-th time step
    STATES = zeros(L, t0)
    snapshot = zeros(L, L)
    # insert Current measurement
    hop_counter = 0
    # perform L*t0 update steps
    for t in 1:(t0*L)
        state, hop_counter = update(state, hop_counter, α, β, p)
        # save snapshot
        if t%L == 0
            STATES[:, div(t, L)] = vec(mean(snapshot, dims=2))
        else
            snapshot[:, t%L] = state 
        end
    end
    return STATES, hop_counter/10_000
end