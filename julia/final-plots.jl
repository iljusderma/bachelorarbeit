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
        label="d=$p2",  
        ylims=[0, 1],
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=:outerright) #aspectratio=500)
    y = zeros(L) .+ α/p2
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{p_2+1}")
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{p_2}{p_2+1}")
    p2 = 0.85
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    scatter!(densityprofile, msw=0, ms=2,  
        label="d=$p2")
    hline!([α/p2], ls=:dash, label=L"\frac{\alpha}{p2}")
    end

    if phase==2
    # MC Phase
    α=0.8
    p2=0.25
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    plot2 = scatter(densityprofile, msw=0, ms=2,  
        label="d=$p2", 
        title="α=$α, β=$β", titleloc=:left, titlefont=12, 
        ylims=[0, 1],
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=:outerright)
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{p_2+1}")
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{p_2}{p_2+1}")
    end
    
    if phase==3
    # HD Phase
    α = 0.8
    β = 0.4
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    plot3 = scatter(densityprofile, msw=0, ms=2,  
        label="d=$p2", 
        title="α=$α, β=$β", titleloc=:left, titlefont=12, 
        ylims=[0, 1],
        ylabel=L"\langle \rho_i \rangle", 
        xlabel="Lattice site i", legend=:outerright)
    y = zeros(L) .+ (1-β/p2)
    plot!(1:251, zeros(251).+ 1/(p2+1), ls=:dash, label=L"\frac{1}{p_2+1}")
    plot!(251:500, zeros(250).+ p2/(p2+1), ls=:dash, label=L"\frac{p_2}{p_2+1}")
    
    p2 = 0.85
    STATES, CURRENT = simulate(α, β, L, t0, p1, p2)
    densityprofile = vec(mean(STATES, dims=2))
    scatter!(densityprofile, msw=0, ms=2,  
        label="d=$p2")
    hline!([1-β/p2], ls=:dash, label=L"\frac{\alpha}{p2}")
    end
    # plot(plot1, plot2, plot3, layout=(1, 3), size=(600,210))
    savefig("XX-densityprofile.pdf")
end

modified_TASEP_phases(3)