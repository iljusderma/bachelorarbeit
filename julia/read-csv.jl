using CSV, Tables, Plots, LaTeXStrings, PlotThemes
theme(:lime)

function calc_current(alpha, beta)
        # contourlines
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

function plot_current(path)
        FLUX = CSV.read(path, Tables.matrix, header=0)
        gridsize = 200
        ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
        heatmap(ALPHA, BETA, FLUX, title= L"Current $J$", 
                xlabel=L"Exit rate $\beta$", 
                ylabel=L"Entry rate $\alpha$")
        Z = @. calc_current(ALPHA, BETA')
        contour!(ALPHA, BETA, Z, levels=[0.05, 0.1, 0.15, 0.2, 0.25], color=:solar, lw=2)
end

function plot_fs_iterations(path)
        DENSITY = CSV.read(path, Tables.matrix, header=0)
        ITERATIONS = range(1000, step=10*1000, length=gridsize)
        scatter(log.(1 ./ ITERATIONS) .*(-1), DENSITY, ms=1, 
                        title="Finite scaling with iterations", 
                        xlabel=L"$-ln \left( \frac{1}{iterations} \right) $", 
                        ylabel=L"\langle \rho \rangle",
                        label="density for iterations ⟶ ∞")
end

DENSITY = CSV.read("/home/ilja/Documents/coding/bachelorarbeit/fs-L-200.csv", Tables.matrix, header=0)
L_array = range(100, step=50, length=gridsize)
scatter(log.(1 ./ L_array) .*(-1), DENSITY, ms=1, 
                        title="Finite scaling with lattice size L", 
                        xlabel=L"$-ln \left( \frac{1}{L} \right) $", 
                        ylabel=L"\langle \rho \rangle",
                        label="density for L ⟶ ∞")