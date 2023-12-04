using CSV, Tables, Plots, LaTeXStrings, PlotThemes
theme(:default)

function expected_current()
        gridsize = 200
        ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
        EXPECTED_CURRENT = @. calc_current(ALPHA, BETA')
        heatmap(ALPHA, BETA, EXPECTED_CURRENT, title= L"Current $J$", 
        xlabel=L"Exit rate $\beta$", 
        ylabel=L"Entry rate $\alpha$")
end

function calc_current(alpha, beta, p2)
        # contourlines
        current = 0
        if alpha <= beta && alpha < 0.5
                current = alpha*(1 - alpha)
        elseif beta < alpha && beta < 0.5
                current = beta*(1 - beta)
        elseif alpha > 0.5 && beta > 0.5
                current = p2/(1+p2)^2
        end
        return current
end

function plot_current_map(path)
        p2 = 0.3
        FLUX = CSV.read(path, Tables.matrix, header=0)
        gridsize = size(FLUX)[1]
        ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
        Z = @. calc_current(ALPHA, BETA', p2)
        h = heatmap(ALPHA, BETA, FLUX, title= L"Current $J$", 
                xlabel=L"Exit rate $\beta$", 
                ylabel=L"Entry rate $\alpha$")
        contour!(ALPHA, BETA, Z, levels=range(0, p2/(1+p2)^2, 6)[2:end], color=:lightrainbow, lw=3)
        hline!([p2/(1+p2)], ls=:dash, color=:black, label=L"\alpha_c", legend=:topright)
        vline!([p2/(1+p2)], ls=:dash, color=:black, label=L"\beta_c")
        display("image/png", h)
end

function plot_current_line(path)
        CURRENT = CSV.read(path, Tables.matrix, header=0)
        gridsize = 200
        index = 51              # alpha fest, beta variabel
        CURRENT_line = CURRENT[:, index]
        ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
        alpha = round(ALPHA[index], digits=2)
        scatter(BETA, CURRENT_line, title="Current for alpha = $alpha", 
                xlabel=L"Exit rate $\beta$", 
                ylabel=L"Current $J$", ms=1,
                legend=false)
end
function plot_alpha_rho(path)
        RHO = CSV.read(path, Tables.matrix, header=0)
        gridsize = 200
        index = 51              # alpha fest, beta variabel
        CURRENT_line = CURRENT[:, index]
        ALPHA, BETA = range(0, 1, gridsize), range(0, 1, gridsize)
        alpha = round(ALPHA[index], digits=2)
        scatter(BETA, CURRENT_line, title="Current for alpha = $alpha", 
                xlabel=L"Exit rate $\beta$", 
                ylabel=L"Current $J$", ms=1,
                legend=false)
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

plot_current_map("/home/ilja/bachelorarbeit/current-200-impurity.csv")