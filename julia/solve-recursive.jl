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
    rho1, rho2 = 0, 0
    alpha2 = calc_alpha2(p2, rho1)
    beta1 = calc_beta1(p2, rho2)
    rho1 = calc_rho(α, beta1)
    rho2 = calc_rho(alpha2, β)
    
    for i in 1:100_000
        alpha2 = calc_alpha2(p2, rho1)
        beta1 = calc_beta1(p2, rho2)
        #=
        println("------")
        println(rho1)
        println(rho2)
        =#
        if rho1 == calc_rho(α, beta1) && rho2 == calc_rho(alpha2, β)
            println("i=$i")
            rho1 = round(calc_rho(α, beta1), digits=4)
            rho2 = round(calc_rho(alpha2, β), digits=4)
            println(rho1)
            println(rho2)
            break
        elseif i == 100_000
            println("Keine Konvergenz")
        else
            rho1 = calc_rho(α, beta1)
            rho2 = calc_rho(alpha2, β)
        end
    end
end

iterate(0.8, 0.8, 0.3)