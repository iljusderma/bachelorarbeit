include("main.jl")
using CSV, Tables, LsqFit

function determine_current_map()
    # counts the number of hops at a site and devide by iterations
    gridsize = 200
    CURRENT = zeros(gridsize, gridsize)
    ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
    for i in 1:length(CURRENT)
        println(i)
        if (i)%10 == 0
            println(round(i/gridsize^2*100), "%")
        end
        # index to row x column
        a, b = (i - 1)%gridsize + 1, div(i-1, gridsize) + 1
        α, β = ALPHA[a], BETA[b]
        STATES, CURRENT[i] = simulate(α, β, L, t0)
        CSV.write("current-"*string(gridsize)*".csv",  Tables.table(CURRENT), writeheader=false)
    end
end

function mean_density(α, β, L, t0, p)
    nos = 3 # Nos: number of simulations
    @time begin
    TOTAL_DENSITY = zeros(Float64, t0, nos)
    for i in 1:nos
        println(i/nos*100, "%")
        STATES, curr = simulate(α, β, L, t0, p)
        TOTAL_DENSITY[:, i] = vec(mean(STATES, dims=1))
    end
    total_density = vec(mean(TOTAL_DENSITY, dims=2))
    end

    # exponential fit
    # @. model(x, p) = p[1]*exp(-p[2]*(x-p[3])) + p[4] -> p[3] ist fast null 
    @. model(x, par) = par[1]*exp(-par[2]*p^2*x) + par[3]
    xdata = (1:t0) ./ t0
    ydata = total_density
    par0 = [0.2, 11.4, 0.65]
    fit = curve_fit(model, xdata, ydata, par0)
    params = @. round(fit.param, digits=2)
    println("Fit-Parameter:", params)
    # determine approached density
    interval = floor(Int, 0.35*t0)# define interval
    mean_total_density = mean(total_density[interval:end])
    println("Approach to density = $mean_total_density")
    # uncertainty of the approach
    error = std(TOTAL_DENSITY[interval:end, :])
    println("with error = $error")
    gr()
    # plot data
    scatter((1:t0) ./ t0, total_density, label="Calculated density evolution", 
            ms=0.1, title="Mean over $nos simulations", 
            xlabel="Standardized time t", 
            ylabel=L"\langle \overline{\rho} \rangle", ylims=(0, 1))
    # plot fit
    plot!(xdata, model(xdata, params), label=L"Fit with $ae^{-bt}+c$")
    annotate!((0.5, 0.4), "(α, β, p) = ($α, $β, $p)")
end

function small_plot()
    # taken by hand with mean_density
    p = [1, 0.95, 0.9, 0.85, 0.8, 0.75]
    density = [0.293, 0.314, 0.331, 0.352, 0.371, 0.395]
    density_err = [0.031, 0.033, 0.033, 0.034, 0.034, 0.034]
    gr()
    scatter(p, density, yerr=density_err, legend=:bottomright, 
        label="Calculated mean density with errorbars", 
        title=L"$\overline{\rho}(p)$-diagram", 
        xlabel="Hop rate p", 
        ylabel=L"Mean density $\overline{\rho}$")
    hline!([0.3], ls=:dash, 
    label="Expected density value")
end

t0 = 10_000 # one time unit includes L updates of the lattice
L = 200

# mean_density(0.3, 0.6, L, t0, 1) # (α, β, L, t0, p)
# determine_current_map()
small_plot()