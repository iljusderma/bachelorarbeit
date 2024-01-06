include("main.jl")
using CSV, Tables, LsqFit

function determine_current_map(L, t0)
    gridsize = 200
    CURRENT = zeros(gridsize, gridsize)
    ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
    for i in 1:length(CURRENT)
        if (i)%10 == 0
            println(round(i/gridsize^2*100), "%")
        end
        # index to row x column
        a, b = (i - 1)%gridsize + 1, div(i-1, gridsize) + 1
        α, β = ALPHA[a], BETA[b]
        STATES, CURR = simulate(α, β, L, t0, 1, 0.3)
        CURRENT[i] = mean(CURR[1_000:end])
        CSV.write("current-"*string(gridsize)*"-impurity.csv",  Tables.table(CURRENT), writeheader=false)
    end
end

function calc_V_p_curve(ALPHA, BETA, gridsize, t0, L)
    # phase transition
    RHO = zeros(gridsize)
    println("start")
    for i in 1:gridsize
        if i%5 ==0
            println(round(i/gridsize*100), "%")
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

function V_p_1order(t0)
    gridsize = 50
    α0 = 0.1
    ALPHA = range(0, 2*α0, gridsize)
    BETA = .- ALPHA .+ 2*α0  # apply shock at alpha = beta = α0
    RHO10 = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 10)
    RHO500 = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 500)
    p = BETA .- ALPHA
    scatter(p, RHO10, label="L=10", 
        xlabel=L"p", ylabel=L"V", ms=3, 
        msw=0, ylims=[0, 1], legend_font=12)
    scatter!(p, RHO500, label="L=500", ms=3, msw=0)
    leap = zeros(gridsize)
    mid = Int(floor(gridsize*0.5))
    leap[1:mid] .= -0.5 .*p[1:mid] .+ α0
    leap[mid+1:end] .= -0.5 .*p[mid+1:end] .+ (1-α0)
    plot!(p[1:25], leap[1:25], lw=2, label="L ⟶ ∞")
    plot!(p[25:26], leap[25:26], linestyle=:dash, lw=2, label=false, color=palette(:default)[3])
    plot!(p[26:50], leap[26:50], lw=2, label=false, color=palette(:default)[3])
    savefig("plot.pdf")
    CSV.write("V-p-1order.csv",  Tables.table([p, RHO10, RHO500]), writeheader=false)
end

function V_p_2order(t0)
    gridsize = 50
    α0 = 0.7
    ALPHA = zeros(gridsize) .+ α0
    BETA = range(0.2, 0.8, gridsize) # apply transition at β=0.5
    RHO10 = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 10)
    RHO500 = calc_V_p_curve(ALPHA, BETA, gridsize, t0, 500)
    p = BETA .- ALPHA
    scatter(p, RHO10, label="L=10", 
        xlabel=L"p", ylabel=L"V", 
        msw=0, ylims=[0, 1], legendfont=12)
    scatter!(p, RHO500, label="L=500", ms=3, msw=0)
    # draw limit L → ∞
    y = zeros(gridsize)
    mid = Int(floor(0.5*length(y)))
    y[1:mid] .= 1 .- (p[1:mid] .+ α0) # ∼1-β
    y[mid+1:end] .= 0.5
    plot!(p, y, label="L ⟶ ∞", lw=2)
    savefig("plot.pdf")
    CSV.write("V-p-2order.csv",  Tables.table([p, RHO10, RHO500]), writeheader=false)
end

function mean_density(α, β, L, t0, p)
    nos = 500 # Nos: number of simulations
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

function rholeft_p2(α, β, L, t0)
    # plot rho_leftbulk over p2 to 
    # find critical p_2 describing a phase transition
    n = 100
    P2 = range(0, 1, n)
    RHOLEFT = zeros(n)
    for (i, p2) in enumerate(P2)
        STATES, CURRENT = simulate(α, β, L, t0, 1, p2)
        cut_STATES = STATES[:, 3_000:end]
        densityprofile = vec(mean(cut_STATES, dims=2))
        RHOLEFT[i] = abs(α - mean(densityprofile[round(Int, L/8):round(Int, 3*L/8)]))
        println(p2)
    end
    DATA = [vec(P2)'; vec(RHOLEFT)']
    CSV.write("rholeft-p2.csv",  Tables.table(DATA), writeheader=false)
    p = scatter(P2, RHOLEFT, 
        xlabel=L"p_2", 
        ylabel=L"\langle \rho_{left} \rangle - α",
        label="α=$α, β=$β, L=$L")
    display("image/png", p) # export as png
end

function rhoright_d(α, β, L, t0)
    # plot rho_leftbulk over p2 to 
    # find critical p_2 describing a phase transition
    n = 100
    D = range(0, 1, n)
    RHOright = zeros(n)
    for (i, d) in enumerate(D)
        STATES, CURRENT = simulate(α, β, L, t0, 1, d)
        cut_STATES = STATES[:, 3_000:end]
        densityprofile = vec(mean(cut_STATES, dims=2))
        RHOright[i] = abs(α - mean(densityprofile[round(Int, 5/8*L):round(Int, 7/8*L)]))
        println(d)
    end
    DATA = [vec(D)'; vec(RHOright)']
    CSV.write("RHOright-d.csv",  Tables.table(DATA), writeheader=false)
    p = scatter(D, RHOright, 
        xlabel=L"d", 
        ylabel=L"\langle \rho_{right} \rangle - α",
        label="α=$α, β=$β, L=$L")
    savefig("plot.pdf")
end

function J_d(α, β, L, t0)
    # find critical p_2 describing a phase transition
    n = 100
    D = range(0, 1, n)
    J = zeros(n)
    for (i, d) in enumerate(D)
        STATES, CURRENT = simulate(α, β, L, t0, 1, D)
        J[i] = mean(CURRENT) - 0.25
        println(p2)
    end
    DATA = [vec(P2)'; vec(J)']
    CSV.write("J-d.csv",  Tables.table(DATA), writeheader=false)
    p = scatter(D, J, 
        xlabel=L"d", 
        ylabel=L"J - \frac{1}{4}",
        label="α=$α, β=$β, L=$L")
    display("image/png", p) # export as png
end

t0 = 20_000 # one time unit includes L updates of the lattice
L = 200

# determine_current_map(L, t0)
# small_plot()
# mean_density(0.2, 0.6, L, t0, 0.9) # (α, β, L, t0, p)
# rholeft_p2(0.3, 0.8, 200, 20_000)
# rhoright_d(0.3, 0.8, 200, 20_000)
# J_d(0.4, 0.8, 500, 50_000)
# V_p_1order(50_000)
V_p_2order(100_000)