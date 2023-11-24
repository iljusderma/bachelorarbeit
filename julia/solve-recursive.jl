using LaTeXStrings

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

function iterate(α, β, p2)
    rhoA, rhoB = 1, 1
    RHOn = zeros((2, 100_000))

    for n in 1:100_000
        alphaB = calc_alpha2(p2, rhoA)
        betaA = calc_beta1(p2, rhoB)
        #=
        println("------")
        println(rho1)
        println(rho2)
        =#
        if rhoA == calc_rho(α, betaA) && rhoB == calc_rho(alphaB, β)
            println("i=$n")
            rhoA = round(calc_rho(α, betaA), digits=4)
            rhoB = round(calc_rho(alphaB, β), digits=4)
            #println(rho1)
            #println(rho2)
            break
        elseif n == 100_000
            println("Keine Konvergenz")
        else
            RHOn[:, n] = vec([rhoA, rhoB])
            rhoA = calc_rho(α, betaA)
            rhoB = calc_rho(alphaB, β)
        end
    end
    RHOAn=RHOn[1, :][findall(x->x!=0, RHOn[1, :])]
    return vec([rhoA, rhoB]), RHOAn
end

a, b = iterate(0.8, 0.8, 0.7)
plot(b[1:25])


#=
α = 0.8
β = 0.8
solution = 100
RHO = zeros((2,solution))
P2 = range(0, 1, solution)
for (i, p2) in enumerate(P2)
    RHO[:, i] = iterate(α, β, p2)
end




# plot rho(p2)-diagram
p = plot(P2, RHO[1,:], 
    xlabel=L"p_2", ylabel=L"\rho", title=L"\alpha=\beta=0.8", 
    label=L"\rho^A", 
    xtickfont=12, ytickfont=12, 
    guidefont=18, legendfont=18,
    size=(800, 500))
plot!(P2, RHO[2, :], label=L"\rho^B")
display("image/png", p)=#