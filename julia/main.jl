using Random, Plots, Statistics, LaTeXStrings

function initialize_state(L)
    # initialize state
    state = bitrand(L)
    return state
end

function update(state, hop_counter, α, β, p1, p2)
    L = length(state)
    site = rand(0:L)
    if site == 0
        if state[1] == 0 && rand() < α
            # injection
            state[1] = 1
        end
    elseif site == L 
        if state[L] == 1 && rand() < β
            # ejection
            state[L] = 0
        end
    elseif site == L/2 + 1
        if state[site] == 1 && state[site + 1] == 0 && rand() < p2
            # hop to right with p2
            state[site + 1] = 1
            state[site] = 0
        end
    else
        if state[site] == 1 && state[site + 1] == 0 && rand() < p1
            # hop to right
            state[site + 1] = 1
            state[site] = 0
            if site == L/2
                hop_counter += 1
            end
        end
    end
    return state, hop_counter
end

function simulate(α, β, L, t0, p1=1, p2=1)
    # initialize state
    state = initialize_state(L)
    # save state history in all_states every n-th time step
    # for current measurement t0/1000
    n = Int(round(t0/t0))
    STATES = BitArray(undef, (L, Int(round(t0/n))))
    # insert Current measurement
    hop_counter = 0
    CURRENT = zeros(Int(round(t0/n)))
    # perform L*t0 update steps
    for t in 1:(t0*L)
        state, hop_counter = update(state, hop_counter, α, β, p1, p2)
        # save snapshot
        if t%(n*L) == 0
            STATES[:, div(t, n*L)] = state
            # determine current every 1_000*L-th update
            CURRENT[div(t, n*L)] = hop_counter/n
            hop_counter = 0
        end
    end
    return STATES, CURRENT
end