using Plots, Statistics, Random, LaTeXStrings

@time begin
function initialState(L::Int64, ratio::Int64)
    state = zeros(L)
    state[1:convert(Int64, ratio*L)] .= + 1
    return shuffle(state)
end

function eject(state::Vector, beta:Int64)
    if state[end] == 1 && rand() < beta
        state[end] = 0
        return state
    end
end

function hop(state::Vector, p::Int64, index::Int64)
    if index%2 == 0
        rand_vals = rand(length(state))
        for i in 1:(length(state)-1)
            if i%2 == 0 && state[i] == 1 && state[i + 1] == 0 && rand_vals[i] < p
                state[i] -= 1
                state[i + 1] += 1
            end
        end
    else
        rand_vals = rand(length(state))
        for i in 1:(length(state)-1)
            if i%2 == 1 && state[i] == 1 && state[i + 1] == 0 && rand_vals[i] < p
                state[i] -= 1
                state[i + 1] += 1
            end
        end
    end
    return state
end

function inject(state, alpha)
    if state[1] == 0 && rand() < alpha
        state[1] = 1
        return state
    end
end

L = 500
iterations = 10*1000
occupied_ratio = 0.5

rates = [0.3, 0.8, 1, 0]
state = initialState(L, occupied_ratio)
all_states = zeros(Float64, iterations, L)

for i in 1:iterations
    eject(state, rates[2])
    hop(state, rates[3], i)
    inject(state, rates[1])
    all_states[i, :] = state
end
println("done")
densityprofile = vec(mean(all_states, dims=1))
display(densityprofile)
scatter(1:L, densityprofile, label="$rates", title="Mean density per lattice site")
ylabel!(L"\langle \rho \rangle_i")
xlabel!("Lattice site "*L"i")
end