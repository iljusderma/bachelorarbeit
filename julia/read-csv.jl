using CSV, Tables, Plots, LaTeXStrings, PlotThemes
theme(:lime)

FLUX = CSV.read("/home/ilja/Documents/coding/bachelorarbeit/current-200.csv", Tables.matrix, header=0)
gridsize = 200
ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
heatmap(ALPHA, BETA, FLUX, title= L"Current $J$", 
        xlabel=L"Exit rate $\beta$", 
        ylabel=L"Entry rate $\alpha$")

# contourlines
function calc_current(alpha, beta)
        current = 0
        if alpha < beta && alpha < 0.5
                current = alpha*(1 - alpha)
        end
        if beta < alpha && beta < 0.5
                current = beta*(1 - beta)
        end
        if alpha > 0.5 && beta > 0.5
                current = 0.25
        end
        return current
end

Z = @. calc_current(ALPHA, BETA')
# contour!(ALPHA, BETA, Z, levels=[0.05, 0.1, 0.15, 0.2, 0.25], color=:solar, lw=2)