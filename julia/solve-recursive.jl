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
    solution = 1000
    RHO = zeros((2,solution))
    P2 = range(0, 1, solution)
    for (i, p2) in enumerate(P2)
        RHOAn, RHOBn = iterate(α, β, p2)
        RHO[:, i] = [RHOAn[end], RHOBn[end]]
    end
    
    # plot rho(p2)-diagram
    p = plot(P2, RHO[1,:], 
    xlabel=L"d", ylabel=L"\rho", title="α=$α, β=$β", 
    label=L"\rho^A", titleloc=:left, legendfont=14)
    plot!(P2, RHO[2, :], label=L"\rho^B")
    return p
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

# RHOAn, RHOBn = iterate(0.4, 0.8, 0.6)
# p = plot(RHOBn[1:10])

p = rho_p2_diagram(0.4, 0.8)
savefig("plot.pdf")