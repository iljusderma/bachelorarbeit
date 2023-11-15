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

function iterate(α, p2, β)
    rho1, rho2 = 0.4, 0.5
    alpha2 = calc_alpha2(p2, rho1)
    beta1 = calc_beta1(p2, rho2)
    rho1 = calc_rho(α, beta1)
    rho2 = calc_rho(alpha2, β)
    
    for i in 1:100
        alpha2 = calc_alpha2(p2, rho1)
        beta1 = calc_beta1(p2, rho2)
        println("------")
        println(rho1)
        println(rho2)
        if rho1 == calc_rho(α, beta1) && rho2 == calc_rho(alpha2, β)
            print("i=$i")
            break
        else
            rho1 = calc_rho(α, beta1)
            rho2 = calc_rho(alpha2, β)
        end
    end
end

iterate(0.3, 0.4, 0.6)