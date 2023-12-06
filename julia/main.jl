using Random, Plots, Statistics, LaTeXStrings

function initialize_state(L)
    # initialize state
    state = zeros(L)
    state[1:Int(floor(L*0.5))] .= 1
    return shuffle(state)
end

function pole_update(state, hop_counter, α, β, p1, p2)
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
    elseif site == L/2 + 1
        if state[site] == 1 && state[site + 1] == 0 && rand() < p2
            # hop to right with p2
            state[site + 1] += 1
            state[site] -= 1
        end
    else
        if state[site] == 1 && state[site + 1] == 0 && rand() < p1
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

function simulate(α, β, L, t0, p1=1, p2=1)
    # initialize state
    state = initialize_state(L)
    # save state history in all_states every n-th time step
    n = L
    STATES = zeros(L, Int((t0*L)/n))
    # insert Current measurement
    hop_counter = 0
    CURRENT = zeros(Int((t0*L)/n))
    # perform L*t0 update steps
    for t in 1:(t0*L)
        state, hop_counter = pole_update(state, hop_counter, α, β, p1, p2)
        # save snapshot
        if t%n == 0
            STATES[:, div(t, n)] = state
            CURRENT[div(t, n)] = hop_counter
            hop_counter = 0
        end
    end
    return STATES, CURRENT
end