include("main.jl")
using CSV, Tables


function fs_iterations(L, occupied_ratio, gridsize, flux_steps, rates)
    DENSITY = zeros(gridsize)
    T = range(1000, step=10*1000, length=gridsize)
    @time begin
    for (index, t) in enumerate(T)
        if (index)%10 == 0
            println(round(index/gridsize*100, digits=2), "%")
        end
        # perform single simulation
        state = initialState(L, occupied_ratio)
        all_states, FLUX = SLUpdate(t, state, rates, flux_steps)
        # determine densityprofile
        cut_states = all_states[:, L:t]
        densityprofile = vec(mean(cut_states, dims=2))
        # determine average density
        DENSITY[index] = mean(densityprofile[1:(L-10)])
        CSV.write("fs-iterations-"*string(gridsize)*".csv",  Tables.table(DENSITY), writeheader=false)
    end
    end
    scatter(1 ./ T, DENSITY)
end

function fs_L(t, occupied_ratio, gridsize, flux_steps, rates)
    DENSITY = zeros(gridsize)
    L_array = range(100, step=100, length=gridsize)
    for (index, L) in enumerate(L_array)
        if (index)%10 == 0
            println(round(index/gridsize*100, digits=2), "%")
        end
        # perform single simulation
        state = initialState(L, occupied_ratio)
        all_states, FLUX = SLUpdate(t, state, rates, flux_steps)
        # determine densityprofile
        cut_states = all_states[:, L:t]
        densityprofile = vec(mean(cut_states, dims=2))
        # determine average density
        DENSITY[index] = mean(densityprofile[1:(L-10)])
        CSV.write("fs-L-"*string(gridsize)*".csv",  Tables.table(DENSITY), writeheader=false)
    end
    scatter(1 ./ L_array, DENSITY)
end

function impurity_MC_density_deviation(α, β, p2)
    # n: number of different simulations
    n = 50
    # DATA: first line includes values of L and second line the deviations
    DATA = zeros(2, n)
    # generate different even L in log range
    DATA[1, :] = round.(range(100, 2000, n)) .* 2
    for (i, L) in enumerate(DATA[1, :])
        STATES, CURRENT = simulate(α, β, Int(L), 50_000, 1, p2)
        cut_STATES = STATES[:, 10_000:end]
        densityprofile = vec(mean(cut_STATES, dims=2))
        # difference between the approximated density and the average density in the middle of a sublattice
        DATA[2, i] = 1/(p2+1) - mean(densityprofile[round(Int, L/8):round(Int, 3*L/8)])
        println(L)
    end
    CSV.write("fs-impurity-MC.csv",  Tables.table(DATA), writeheader=false)
    scatter(DATA[1,:], DATA[2,:])
end

function impurity_MC_current_deviation(α, β, p2)
    # n: number of different simulations
    n = 6
    # DATA: first line includes values of L and second line the deviations
    DATA = zeros(3, n)
    # generate different even L in log range
    DATA[1, :] = round.(2 .^ range(7, 7 + n))
    for (i, L) in enumerate(DATA[1, :])
        STATES, CURRENT = simulate(α, β, Int(L), 1_000_000, 1, p2)
        # difference between the approximated density and the average density in the middle of a sublattice
        DATA[2, i] = mean(CURRENT[3:end]) - p2/(p2+1)^2
        DATA[3, i] = std(CURRENT[3:end])/sqrt(length(CURRENT[3:end]))
        println(L)
    end
    CSV.write("fs-impurity-MC-current-deviation.csv",  Tables.table(DATA), writeheader=false)
    scatter(DATA[1,:], DATA[2,:])
end

impurity_MC_current_deviation(1, 1, 1)