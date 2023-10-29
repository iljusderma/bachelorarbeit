include("main.jl")
using CSV, Tables
using LsqFit

function iterate_fakeflux(rates, flux_steps, t, L, occupied_ratio)
    # counts the number of hops at a site and devide by iterations
    gridsize = 200
    FLUX_Matrix = zeros(Float64, gridsize, gridsize)
    ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
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
            all_states, FLUX = RANDUpdate(t, state, rates, flux_steps)
            popfirst!(FLUX)
            FLUX_Matrix[index] = mean(FLUX)
            CSV.write("flux-"*string(gridsize)*".csv",  Tables.table(FLUX_Matrix), writeheader=false)
        end
    end
    end
    return FLUX_Matrix
end

function iterate_current(CURRENT, rates, flux_steps, t, L, occupied_ratio)
    gridsize = 200
    FLUX_Matrix = zeros(Float64, gridsize, gridsize)
    ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
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
            all_states, FLUX = SLUpdate(t, state, rates, flux_steps)
            # determine J = rho_i(1-rho_ii)
            cut_states = all_states[L:t, :]
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

function mean_density(L, occupied_ratio, flux_steps, rates, t)
    nos = 100 # Nos: number of simulations
    @time begin
    TOTAL_DENSITY = zeros(Float64, t, nos)
    for i in 1:nos
        state = initialState(L, occupied_ratio)
        all_states, CURRENT = SLUpdate(t, state, rates, flux_steps)
        TOTAL_DENSITY[:, i] = vec(mean(all_states, dims=1))
    end
    total_density = vec(mean(TOTAL_DENSITY, dims=2))
    end
    p = scatter(total_density, legend=false, ms=0.1, title="nos = $nos")
    return total_density, p
end

L = 500
t = 10*1000
occupied_ratio = 0.5
flux_steps = 50
rates = [0.6, 0.3, 0.9, 0]
gr()
total_density, p = mean_density(L, occupied_ratio, flux_steps, rates, t)

@. model(x, p) = p[1]*exp(p[2]*(x-p[3])) + p[4]
xdata = 1.0:t
ydata = total_density
p0 = [1.0,1.0,1.0,1.0]

# Überprüfe und bereinige die Daten
inds = .!isnan.(ydata) .& .!isinf.(ydata)  # Indizes ohne NaN oder Inf
xdata_clean = xdata[inds]
ydata_clean = ydata[inds]

fit = curve_fit(model, xdata, ydata, p0)
params = coef(fit)
display(p)
plot!(xdata, model(xdata, params))