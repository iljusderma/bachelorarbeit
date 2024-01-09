using LaTeXStrings, Plots

function calc_alpha2(p2, rho1)
    return p2*rho1    
end

function calc_beta1(p2, rho2)
    return p2*(1-rho2)
end

function calc_rho(α, β)
    if α<β && α<0.5 #LD
        ρ = α
    elseif α>=β && β<0.5 #HD
        ρ = 1 - β
    else
        ρ = 1/2
    end
    return ρ
end

function rho_p2_diagram(α, β)
    # plot theory
    solution = 1000
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    plot(P2, RHO[1,:], lw=2, label="Fixed points")
    plot!(P2, RHO[2,:], lw=2, color=palette(:default)[1], label=false)
    # scatter numerics
    solution = 50
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    scatter!(P2, RHO[1,:], 
        xlabel=L"d", ylabel=L"\rho", title="α=$α, β=$β", 
        label=L"$\rho^A$", titleloc=:left,
        legendfontsize=10, msw=0, ms=4)
    scatter!(P2, RHO[2, :], label=L"$\rho^B$", msw=0, ms=4)
    savefig("plot.pdf")
end

function iterate(α, β, p2)
    rhoA, rhoB = 1, 1
    # save all values in an array (5 is no possible value -> filter later)
    RHOn = zeros((2, 10_000)).+5

    for n in 1:10_000
        alphaB = calc_alpha2(p2, rhoA)
        betaA = calc_beta1(p2, rhoB)
        rhoA_old, rhoB_old = rhoA, rhoB
        rhoA, rhoB = calc_rho(α, betaA), calc_rho(alphaB, β)
        RHOn[:, n] = vec([rhoA_old, rhoB_old])
        rhoA = calc_rho(α, betaA)
        rhoB = calc_rho(alphaB, β)
    end
    RHOAn = RHOn[1, :][findall(x->x!=5, RHOn[1, :])]
    RHOBn = RHOn[2, :][findall(x->x!=5, RHOn[2, :])]
    return RHOAn, RHOBn
end

function final()
    # MC
    α, β = 0.8, 0.8
    # plot theory
    solution = 1000
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    plot1 = plot(P2, RHO[1,:], lw=2, label="Fixed points")
    plot!(P2, RHO[2,:], lw=2, color=palette(:default)[1], label=false)
    # scatter numerics
    solution = 50
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    scatter!(P2, RHO[1,:], 
        xlabel=L"d", ylabel=L"\rho", title="a) α=$α, β=$β", 
        label=L"$\rho^A$", titleloc=:left,
        legendfontsize=10, msw=0, ms=3)
    scatter!(P2, RHO[2, :], label=L"$\rho^B$", msw=0, ms=3, 
        left_margin=2mm)
    # LD 
    α, β = 0.4, 0.8
    # plot theory
    solution = 1000
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    plot2 = plot(P2, RHO[1,:], lw=2, label="Fixed points")
    plot!(P2, RHO[2,:], lw=2, color=palette(:default)[1], label=false)
    # scatter numerics
    solution = 50
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    scatter!(P2, RHO[1,:], 
        xlabel=L"d", ylabel=L"\rho", title="b) α=$α, β=$β", 
        label=L"$\rho^A$", titleloc=:left,
        legendfontsize=10, msw=0, ms=3)
    scatter!(P2, RHO[2, :], label=L"$\rho^B$", msw=0, ms=3, 
        bottom_margin=5mm, left_margin=5mm)
    plot(plot1, plot2, layout=2, size=(800, 400))
    savefig("plot.pdf")
end

# RHOAn, RHOBn = iterate(0.4, 0.8, 0.6)
# p = plot(RHOBn[1:10])

# rho_p2_diagram(0.4, 0.8)
final()