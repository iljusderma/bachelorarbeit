using Plots

gridsize = 200
rates = [0.6, 0.3, 1, 0]
FLUX_Matrix = zeros(Float64, gridsize, gridsize)
ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)

for (index, flux) in enumerate(FLUX_Matrix)
    if flux == 0
        if (index)%10 == 0
            println(round(index/gridsize^2*100), "%")
        end
        # index to row x column
        a, b = (index - 1)%gridsize + 1, div(index-1, gridsize) + 1
        rates[1], rates[2] = ALPHA[a], BETA[b]
        if rates[1]<rates[2] && rates[1] < 0.5
            FLUX_Matrix[index] = rates[1]*(1-rates[1])
        end
        if rates[2]<rates[1] && rates[2] < 0.5
            FLUX_Matrix[index] = rates[2]*(1-rates[2])
        end
        if rates[1]>0.5 && rates[2]>0.5
            FLUX_Matrix[index] = 0.25
        end
    end
end
display(FLUX_Matrix)
heatmap(ALPHA, BETA, FLUX_Matrix, title= L"Current $J$", 
    xlabel=L"Exit rate $\beta$", 
    ylabel=L"Entry rate $\alpha$")