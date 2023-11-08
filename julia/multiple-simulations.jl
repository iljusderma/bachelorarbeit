include("main.jl")
using CSV, Tables, LsqFit

function determine_current_map()
    # counts the number of hops at a site and devide by iterations
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

t0 = 10_000 # one time unit includes L updates of the lattice
L = 200

determine_current_map()