include("main.jl")
using CSV, Tables


function fs_iterations(L, occupied_ratio, gridsize, flux_steps, rates)
    DENSITY = zeros(gridsize)
    ITERATIONS = range(1000, step=10*1000, length=gridsize)
    @time begin
    for (index, iterations) in enumerate(ITERATIONS)
        if (index)%10 == 0
            println(round(index/gridsize*100, digits=2), "%")
        end
        # perform single simulation
        state = initialState(L, occupied_ratio)
        all_states, FLUX = SLUpdate(iterations, state, rates, flux_steps)
        # determine densityprofile
        cut_states = all_states[L:iterations, :]
        densityprofile = vec(mean(cut_states, dims=1))
        # determine average density
        DENSITY[index] = mean(densityprofile[1:(L-10)])
        CSV.write("fs-iterations-"*string(gridsize)*".csv",  Tables.table(DENSITY), writeheader=false)
    end
    end
    scatter(1 ./ ITERATIONS, DENSITY)
end

function fs_L(iterations, occupied_ratio, gridsize, flux_steps, rates)
    DENSITY = zeros(gridsize)
    L_array = range(100, step=100, length=gridsize)
    for (index, L) in enumerate(L_array)
        if (index)%10 == 0
            println(round(index/gridsize*100, digits=2), "%")
        end
        # perform single simulation
        state = initialState(L, occupied_ratio)
        all_states, FLUX = SLUpdate(iterations, state, rates, flux_steps)
        # determine densityprofile
        cut_states = all_states[L:iterations, :]
        densityprofile = vec(mean(cut_states, dims=1))
        # determine average density
        DENSITY[index] = mean(densityprofile[1:(L-10)])
        CSV.write("fs-L-"*string(gridsize)*".csv",  Tables.table(DENSITY), writeheader=false)
    end
    scatter(1 ./ L_array, DENSITY)
end

stat_L = 100
stat_iterations = 50*1000
occupied_ratio = 0.5
gridsize = 200
flux_steps = 50
rates = [0.3, 0.8, 1, 0]

# fs_iterations(stat_L, occupied_ratio, gridsize, flux_steps, rates)
fs_L(stat_iterations, occupied_ratio, gridsize, flux_steps, rates)