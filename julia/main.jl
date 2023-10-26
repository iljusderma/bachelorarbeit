using Plots, Statistics, Random, LaTeXStrings
using PlotThemes
theme(:juno)

function initialState(L::Int64, ratio::Float64)
    state = zeros(L)
    state[1:convert(Int64, ratio*L)] .= + 1
    return shuffle(state)
end

function eject(state::Vector, beta::Float64)
    if state[end] == 1 && rand() < beta
        state[end] = 0
        return state
    end
end

function SLhop(state::Vector, p::Float64, index::Int64, hop_counter)
	L = length(state)
    if index%2 == 0
        rand_vals = rand(L)
        for i in 1:(L-1)
            if i%2 == 0 && state[i] == 1 && state[i + 1] == 0 && rand_vals[i] < p
                state[i] -= 1
                state[i + 1] += 1
				if i == L/2
					hop_counter += 1
				end
            end
        end
    else
        rand_vals = rand(L)
        for i in 1:(L-1)
            if i%2 == 1 && state[i] == 1 && state[i + 1] == 0 && rand_vals[i] < p
                state[i] -= 1
                state[i + 1] += 1
				if i == L/2
					hop_counter += 1
				end
            end
        end
    end
    return state, hop_counter
end

function RANDhop(state::Vector, p::Float64, index::Int64, hop_counter)
	L = length(state)
	# sequence in random order
	sequence = shuffle(1:(L-1))
	# generate random values for p
	rand_vals = rand(L)
	# loop through all sites in random order (sequence)
	for site in sequence
		if state[site] == 1 && state[site + 1] == 0 && rand_vals[site] < p
			state[site] -= 1
            state[site + 1] += 1
			if site == div(L, 2)
				hop_counter += 1
			end
		end
	end
    return state, hop_counter
end

function inject(state, alpha)
    if state[1] == 0 && rand() < alpha
        state[1] = 1
        return state
    end
end

function SLUpdate(iterations, state, rates, flux_steps)
	all_states = zeros(Float64, iterations, length(state))
	FLUX = zeros(div(iterations, flux_steps))
	hop_counter = 0
	for i in 1:iterations
		eject(state, rates[2])
		state, hop_counter = SLhop(state, rates[3], i, hop_counter)
		inject(state, rates[1])
		all_states[i, :] = state
		if i%flux_steps == 0
			FLUX[div(i, flux_steps)] = hop_counter/flux_steps
			hop_counter = 0
		end
	end
	return all_states, FLUX
end

function RANDUpdate(iterations, state, rates, flux_steps)
	all_states = zeros(Float64, iterations, L)
	FLUX = zeros(div(iterations, flux_steps))
	hop_counter = 0
	for i in 1:iterations
		eject(state, rates[2])
		state, hop_counter = RANDhop(state, rates[3], i, hop_counter)
		inject(state, rates[1])
		all_states[i, :] = state
		if i%flux_steps == 0
			FLUX[div(i, flux_steps)] = hop_counter/flux_steps
			hop_counter = 0
		end
	end
	return all_states, FLUX
end