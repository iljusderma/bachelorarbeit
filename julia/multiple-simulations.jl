include("main.jl")
using CSV, Tables, LaTeXStrings
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
    nos = 500 # Nos: number of simulations
    @time begin
    TOTAL_DENSITY = zeros(Float64, t, nos)
    for i in 1:nos
        state = initialState(L, occupied_ratio)
        all_states, CURRENT = SLUpdate(t, state, rates, flux_steps)
        TOTAL_DENSITY[:, i] = vec(mean(all_states, dims=1))
    end
    total_density = vec(mean(TOTAL_DENSITY, dims=2))
    end
    # exponential fit
    # @. model(x, p) = p[1]*exp(-p[2]*(x-p[3])) + p[4] -> p[3] ist fast null 
    @. model(x, p) = p[1]*exp(-p[2]*rates[3]^2*x) + p[3]
    xdata = (1:t) ./ t
    ydata = total_density
    p0 = [0.2, 11.4, 0.65]
    fit = curve_fit(model, xdata, ydata, p0)
    params = fit.param
    println(params[2])
    # determine approached density
    interval = floor(Int, 0.35*t)# define interval
    mean_total_density = mean(total_density[interval:end])
    println("Approach to density = $mean_total_density")
    # uncertainty of the approach
    error = std(TOTAL_DENSITY[interval:end, :])
    println("with error = $error")
    plotlyjs()
    # plot data
    scatter((1:t) ./ t, total_density, label="Calculated density evolution", 
            ms=0.1, title="Mean over $nos simulations, $rates", 
            legend = :outertopright)
    # plot fit
    plot!(xdata, model(xdata, params), label=L"Fit with $ae^{-bt}+c$")
end

function small_plot()
    # taken by hand with mean_density
    p = [1, 0.95, 0.9, 0.85, 0.8, 0.75]
    density = [0.7453, 0.7051, 0.6609, 0.6124, 0.5592, 0.5247]
    density_err = [0.0252, 0.0272, 0.0294, 0.0310, 0.0307, 0.0277]
    gr()
    scatter(p, density, yerr=density_err, legend=:bottomright, 
        label="Calculated mean density with errorbars", 
        title=L"$\overline{\rho}(p)$-diagram", 
        xlabel="Hop rate p", 
        ylabel=L"Mean density $\overline{\rho}$")
    plot!(p, zeros(6) .+ 0.7, ls=:dash, 
        label="Expected density value")
end

L = 200
t = 20*1000
occupied_ratio = 0.5
flux_steps = 50
rates = [0.6, 0.3, 0.75, 0]
# mean_density(L, occupied_ratio, flux_steps, rates, t)
small_plot()