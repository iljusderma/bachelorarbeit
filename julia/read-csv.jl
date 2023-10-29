using CSV, Tables, Plots, LaTeXStrings, PlotThemes
theme(:lime)

function expected_current()
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
end

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
        Z = @. calc_current(ALPHA, BETA')
        heatmap(ALPHA, BETA, FLUX, title= L"Current $J$", 
                xlabel=L"Exit rate $\beta$", 
                ylabel=L"Entry rate $\alpha$")
        #contour!(ALPHA, BETA, Z, levels=[0.05, 0.1, 0.15, 0.2, 0.25], color=:solar, lw=2)
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

function plot_fs_density()
        DENSITY = CSV.read("/home/ilja/Documents/coding/bachelorarbeit/fs-L-200.csv", Tables.matrix, header=0)
        L_array = range(100, step=50, length=gridsize)
        scatter(log.(1 ./ L_array) .*(-1), DENSITY, ms=1, 
                        title="Finite scaling with lattice size L", 
                        xlabel=L"$-ln \left( \frac{1}{L} \right) $", 
                        ylabel=L"\langle \rho \rangle",
                        label="density for L ⟶ ∞")
end

plot_current("/home/ilja/Documents/coding/bachelorarbeit/current-200.csv")