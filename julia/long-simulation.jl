include("main.jl")
using CSV, Tables

function iterate_fakeflux(FLUX_Matrix, gridsize, rates, ALPHA, BETA, flux_steps, iterations, L, occupied_ratio)
    # counts the number of hops at a site and devide by iterations
    @time begin
    for (index, flux) in enumerate(FLUX_Matrix)
        if flux == 0
            if (index)%10 == 0
                println(round(index/gridsize^2*100), "%")
            end
            # index to row x column
            a, b = (index - 1)%gridsize + 1, div(index-1, gridsize) + 1
            rates[1], rates[2] = ALPHA[a], BETA[b]
            state = initialState(L, occupied_ratio)
            all_states, FLUX = RANDUpdate(iterations, state, rates, flux_steps)
            popfirst!(FLUX)
            FLUX_Matrix[index] = mean(FLUX)
            CSV.write("flux-"*string(gridsize)*".csv",  Tables.table(FLUX_Matrix), writeheader=false)
        end
    end
    end
    return FLUX_Matrix
end

function iterate_current(CURRENT, gridsize, rates, ALPHA, BETA, flux_steps, iterations, L, occupied_ratio)
    @time begin
    for (index, flux) in enumerate(FLUX_Matrix)
        if flux == 0
            if (index)%10 == 0
                println(round(index/gridsize^2*100, digits=2), "%")
            end
            # index to row x column
            a, b = (index - 1)%gridsize + 1, div(index-1, gridsize) + 1
            # set new alpha, beta
            rates[1], rates[2] = ALPHA[a], BETA[b]
            # perform single simulation
            state = initialState(L, occupied_ratio)
            all_states, FLUX = SLUpdate(iterations, state, rates, flux_steps)
            # determine J = rho_i(1-rho_ii)
            cut_states = all_states[L:iterations, :]
            densityprofile = vec(mean(cut_states, dims=1))
            rho_i = densityprofile[1:L-1] # density at site i
            rho_ii = densityprofile[2:L]  # density at site i + 1
            J = rho_i .* (1 .- rho_ii)
            CURRENT[index] = mean(J)
            CSV.write("current-"*string(gridsize)*".csv",  Tables.table(CURRENT), writeheader=false)
        end
    end
    end
    return CURRENT
end

L = 500
iterations = 5*1000
occupied_ratio = 0.5
gridsize = 200
flux_steps = 50
rates = [0.6, 0.3, 1, 0]
FLUX_Matrix = zeros(Float64, gridsize, gridsize)
ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)

# first intention
# FLUX_Matrix = iterate_fakeflux(FLUX_Matrix, gridsize, rates, ALPHA, BETA, flux_steps, iterations, L, occupied_ratio)

# calculated by Blythe and Evans
FLUX_Matrix =  iterate_current(FLUX_Matrix, gridsize, rates, ALPHA, BETA, flux_steps, iterations, L, occupied_ratio)
display(FLUX_Matrix)
heatmap(FLUX_Matrix)