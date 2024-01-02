include("main.jl")
using PlotThemes, CSV, Tables

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
        xlabel=L"Lattice site $i$", legend=:outerright) #aspectratio=500)
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
        xlabel=L"Lattice site $i$", legend=:outerright)
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
        xlabel=L"Lattice site $i$", legend=:outerright)
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
    h = heatmap(ALPHA, BETA, FLUX, title= L"Current $J$", 
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
    savefig("plot.pdf")
end

# modified_TASEP_phases(3)
plot_current_map_standard("current-200.csv")