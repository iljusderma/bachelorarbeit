include("main.jl")
using CSV, Tables, LsqFit

function determine_current_map()
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

function calc_V_p_curve(ALPHA, BETA, gridsize, t0, L)
    RHO = zeros(gridsize)
    for i in 1:gridsize
        if (i)%10 == 0
            println(round(i/gridsize^2*100), "%")
        end
        # index to row x column
        α = ALPHA[i]
        β = BETA[i]
        STATES, curr = simulate(α, β, L, t0)
        cut_STATES = STATES[:, 4_000:end]
        RHO[i] = mean(cut_STATES)
    end
    return RHO
end

function multipleV_p(t0)
    gridsize = 50
    α0 = 0.1
    ALPHA = range(0, 2*α0, gridsize)
    BETA = .- ALPHA .+ 2*α0  # apply shock at alpha = beta = 0.3
    RHO10 = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 10)
    RHO200 = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 200)
    p = BETA .- ALPHA
    scatter(BETA .- ALPHA, RHO10, label="L=10", 
    title="V(p)-diagram with p(α=β=$α0)=0", xlabel="p", ylabel="V", ms=3, 
    msw=0, ylims=[0, 1])
    scatter!(BETA .- ALPHA, RHO200, label="L=200", ms=3, msw=0)
    leap = zeros(gridsize)
    mid = Int(floor(gridsize*0.5))
    leap[1:mid] .= -0.5 .*p[1:mid] .+ α0
    leap[mid+1:end] .= -0.5 .*p[mid+1:end] .+ (1-α0)
    plot!(p, leap, label="L ⟶ ∞")
    #CSV.write("alpha-rho-"*string(gridsize)*".csv",  Tables.table(RHO), writeheader=false)
end

function V_p(t0)
    gridsize = 50
    trans = 0.1
    ALPHA = range(0, 2*trans, gridsize)
    BETA = .- ALPHA .+ 2*trans  # apply shock at alpha = beta = 0.3
    RHO = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 200)
    p = BETA .- ALPHA
    scatter(p, RHO)
    mid = Int(floor(gridsize*0.5))
    leap[1:mid] .= -0.5 .*p[1:mid] .+ trans
    leap[mid+1:end] .= -0.5 .*p[mid+1:end] .+ 1-trans
    plot!(p, leap, label="L ⟶ ∞")
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

#t0 = 20_000 # one time unit includes L updates of the lattice
#L = 200

# mean_density(0.3, 0.6, L, t0, 1) # (α, β, L, t0, p)
# determine_current_map()
# small_plot()
multipleV_p(20_000)