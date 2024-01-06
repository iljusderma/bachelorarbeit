include("main.jl")
using PlotThemes, CSV, Tables, Measures

function TASEP_phases()
    # initialize lattice parameters
    # lattice size L, injection rate α, ejection rate β, hop rate p
    t0 = 100_000 # one time unit includes L updates of the lattice
    L = 500
    α = 0.4
    β = 0.8
    p1 = 1
    p2 = 1

    # perform update
    @time begin
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    end

    # plot data
    cut_STATES = STATES
    densityprofile = vec(mean(cut_STATES, dims=2))
    totaldensity = vec(mean(STATES, dims=1))

    plot1 = scatter(densityprofile, msw=0, ms=2,  
        label="α=$α, β=$β", title="a) α=$α, β=$β", titleloc=:left, 
        ylims=[0, 1], titlefont=12, 
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=false)

    α=0.8
    @time begin
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    end
    # plot data
    cut_STATES = STATES
    densityprofile = vec(mean(cut_STATES, dims=2))
    totaldensity = vec(mean(STATES, dims=1))

    plot2 = scatter(densityprofile, msw=0, ms=2,  
        label="α=$α, β=$β", title="b) α=$α, β=$β", titleloc=:left, 
        ylims=[0, 1], titlefont=12, 
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=false)

    β=0.4
    @time begin
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    end
    # plot data
    cut_STATES = STATES
    densityprofile = vec(mean(cut_STATES, dims=2))
    totaldensity = vec(mean(STATES, dims=1))

    plot3 = scatter(densityprofile, msw=0, ms=2,  
        label="α=$α, β=$β", title="c) α=$α, β=$β", titleloc=:left, 
        ylims=[0, 1], titlefont=12, 
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=false)

    α=0.3
    β=0.3
    STATES, CURRENT = simulate(α, β, L, 5_000, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))

    plot4 = scatter(densityprofile, msw=0, ms=2,  
        label="α=$α, β=$β", title="d) α=$α, β=$β", titleloc=:left, 
        ylims=[0, 1], titlefont=12, 
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=false)


    #animate(500, t0, STATES, α, β, p1)

    #plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{p_2+1}")
    #plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{p_2}{p_2+1}")
    #println(mean(CURRENT))
    plot(plot1, plot2, plot3, plot4, layout = 4)
    savefig("plot.pdf")
end

function modified_TASEP_phases(phase)
    # initialize lattice parameters
    # lattice size L, injection rate α, ejection rate β, hop rate p
    t0 = Int(1e5) # one time unit includes L updates of the lattice
    L = 500
    α = 0.4
    β = 0.8
    p1 = 1
    p2 = 0.25

    if phase==1
    # LD phase
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    plot1 = scatter(densityprofile, msw=0, ms=2, 
        title="α=$α, β=$β", titleloc=:left, titlefont=12, 
        label=L"$d$ = $0.25$",  
        ylims=[0, 1],
        ylabel=L"\langle \rho_i \rangle", 
        xlabel=L"Lattice site $i$", legend=:outerright, legendfont=12)
    y = zeros(L) .+ α/p2
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{d+1}")
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{d}{d+1}")
    p2 = 0.85
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    scatter!(densityprofile, msw=0, ms=2,  
        label=L"$d$ = $0.85$")
    hline!([α/p2], ls=:dash, label=L"\frac{\alpha}{d}")
    end

    if phase==2
    # MC Phase
    α=0.8
    p2=0.25
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    plot2 = scatter(densityprofile, msw=0, ms=2,  
        label=L"$d$ = $0.25$", 
        title="α=$α, β=$β", titleloc=:left, titlefont=12, 
        ylims=[0, 1],
        ylabel=L"\langle \rho_i \rangle", 
        xlabel=L"Lattice site $i$", legend=:outerright, legendfont=12)
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{d+1}")
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{d}{d+1}")
    p2=0.85
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    scatter!(densityprofile, msw=0, ms=2,  
        label=L"$d$ = $0.85$")
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=false, color=palette(:default)[2])
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=false, color=palette(:default)[3])
    end
    
    if phase==3
    # HD Phase
    α = 0.8
    β = 0.4
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    plot3 = scatter(densityprofile, msw=0, ms=2,  
        label=L"$d$ = $0.25$", 
        title="α=$α, β=$β", titleloc=:left, titlefont=12, 
        ylims=[0, 1],
        ylabel=L"\langle \rho_i \rangle", 
        xlabel=L"Lattice site $i$", legend=:outerright, legendfont=12)
    y = zeros(L) .+ (1-β/p2)
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{d+1}")
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{d}{d+1}")
    p2 = 0.85
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    scatter!(densityprofile, msw=0, ms=2,  
        label=L"$d$ = $0.85$")
    hline!([1-β/p2], ls=:dash, label=L"1-\frac{\beta}{d}")
    end
    # plot(plot1, plot2, plot3, layout=(1, 3), size=(600,210))
    savefig("XX-densityprofile.pdf")
end

function calc_current(alpha, beta, p2)
    # contourlines
    current = 0
    if alpha <= beta && alpha < 0.5
            current = alpha*(1 - alpha)
    elseif beta < alpha && beta < 0.5
            current = beta*(1 - beta)
    elseif alpha > 0.5 && beta > 0.5
            current = p2/(1+p2)^2
    end
    return current
end

function plot_current_map_standard(path)
    FLUX = CSV.read(path, Tables.matrix, header=0)
    gridsize = size(FLUX)[1]
    ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
    Z = @. calc_current(ALPHA, BETA', 1)
    h = heatmap(ALPHA, BETA, FLUX, 
        topmargin=7mm, 
        xlabel=L"Exit rate $\beta$", 
        ylabel=L"Entry rate $\alpha$", 
        aspectratio=true, xlims=(0,1), 
        legendfont=14)
    plot!(range(0.5, 1.0, 10), zeros(10) .+ 0.5, color=:darkgreen, label=false)
    plot!(zeros(10) .+ 0.5, range(0.5, 1.0, 10), color=:darkgreen, label=false)
    plot!(range(0, 0.5, 10),range(0, 0.5, 10), color=:darkgreen, label=false)
    annotate!(0.75, 0.75, ("Maximum current phase", 8, :green, "Helvetica Bold"))
    annotate!(0.25, 0.65, ("High density phase", 8, :green, "Helvetica Bold"))
    annotate!(0.75, 0.25, ("Low density phase", 8, :green, "Helvetica Bold"))
    annotate!(1.1, 1.1, (L"Current $J$", 12))
    savefig("plot.pdf")
end

function plot_critical_d_fromrholeft(path)
    DATA = CSV.read(path, Tables.matrix, header=0)
    A = 0.0:0.05:0.5
    P2 = DATA[1, :]
    P2C = zeros(2, length(A))
    for (i, LEFT) in enumerate(eachrow(DATA[2:end, :]))
            # heap from the right
            P2C[1, i] = P2[(LEFT .- 0.03) .< 0][1]
            # heap from the left
            P2C[2, i] = P2[(LEFT .- 1 ./ (1 .+ LEFT) .+ 0.2 .+ 0.03) .> 0][end]
    end
    xdata = A
    ydata = vec(mean(P2C, dims=1))
    yerr = abs.((P2C[1, :] .- P2C[2, :]) ./ 2)
    p = scatter(xdata, ydata, yerr=yerr, 
            xlabel=L"\alpha", 
            ylabel=L"d_c",
            label="β=0.8, L=200", 
            legend=:outerright)
    # linear fit
    @. model(x, par) = par[1]*x
    par0 = [1.0]
    fit = curve_fit(model, xdata[1:7], ydata[1:7], par0)
    params = @. round(fit.param, digits=2)

    # plot fit
    plot!(xdata[1:7], model(xdata[1:7], params), label=L"Fit with $ax$")
    plot!(xdata[7:end], model(xdata[7:end], params), ls=:dash, label=:none)
    α = 0:0.005:0.5
    plot!(α, α ./ (1 .- α), label=L"\frac{\alpha}{1-\alpha}")

    # first subplot
    # histogram!(
    #     randn(1000),
    #     inset = (1, bbox(0.05, 0.05, 0.25, 0.25, :top, :left)),
    #     ticks = nothing,
    #     subplot = 2,
    #     bg_inside = nothing)
    path = "rholeft-d-02.csv"
    DATA = CSV.read(path, Tables.matrix, header=0)
    scatter!(DATA[1, :], (DATA[2, :]), msw=0, ms=2,
            annotation=(0.7, 0.6, ("α=0.2", 10)), 
            legend=false, 
            inset = (1, bbox(0.02, 0.02, 0.2, 0.3, :top, :left)),
            ticks = nothing,
            subplot = 2,
            bg_inside = nothing)
    path = "rholeft-d-03.csv"
    DATA = CSV.read(path, Tables.matrix, header=0)
    scatter!(DATA[1, :], (DATA[2, :]), msw=0, ms=2,
            annotation=(0.7, 0.53, ("α=0.3", 10)), 
            legend=false, 
            inset = (1, bbox(0.26, 0.02, 0.2, 0.3, :top)),
            ticks = nothing,
            subplot = 3,
            bg_inside = nothing)
    path = "rholeft-d-04.csv"
    DATA = CSV.read(path, Tables.matrix, header=0)
    scatter!(DATA[1, :], (DATA[2, :]), msw=0, ms=2,
            annotation=(0.7, 0.47, ("α=0.4", 10)), 
            legend=false, 
            inset = (1, bbox(0.5, 0.02, 0.2, 0.3, :top)),
            ticks = nothing,
            subplot = 4,
            bg_inside = nothing)

    savefig("plot.pdf")
end

function plot_STATESMAP()
    t0 = Int(5*1e3) # one time unit includes L updates of the lattice
    L = 300
    α = 0.3
    β = 0.3
    p1 = 1
    p2 = 1

    # LD phase
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    heatmap(STATES', c=:grays, 
        title="α=$α, β=$α", titleloc=:left, 
        xlabel=L"Lattice site $i$",
        ylabel=L"Time $t$", legend=false, 
        cbar=false, rightmargin=5mm, 
        legend_font=14)
    # x = [0,  , NaN, 0.5, 1.0]
    # y = [1, 0.4,  NaN, 1.5, 0.8]
    # GR.setarrowsize(1)
    GR.setarrowsize(2)
    plot!([305,50], [4600,2500], lw=3, arrow = (:closed, 2.0), label="(a)")
    annotate!(315,4700, "(a)")
    plot!([305,250], [4200,3000], lw=3, arrow = (:closed, 2.0), label="(b)")
    annotate!(315,4300, "(b)")
    plot!([305,150], [3800,200], lw=3, arrow = (:closed, 2.0), label="(c)")
    annotate!(315,3900, "(c)")
    savefig("plot.pdf")
end

# modified_TASEP_phases(1)
plot_current_map_standard("current-200.csv")
# plot_critical_d_fromrholeft("multiple-rholeft-d.csv")
# plot_STATESMAP()