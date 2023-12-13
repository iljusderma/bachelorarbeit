using CSV, Tables, Plots, LaTeXStrings, LsqFit, Statistics

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
                ylabel=L"Entry rate $\alpha$", 
                aspectratio=true, xlims=(0,1), dpi=500)
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

function plot_fs_impurity_MC_density_deviation(path)
        DATA = CSV.read(path, Tables.matrix, header=0)
        p = scatter((DATA[1, :]) , DATA[2, :],
                title="Deviation in MC phase L→∞", xlabel=L"L",
                ylabel=L"\rho_{left, approx} - \langle \rho_{left} \rangle", 
                label="α=0.8, β=0.8, p2=0.3", dpi=300)
        display("image/png", p)
end

function plot_fs_impurity_MC_current_deviation(path)
        DATA = CSV.read(path, Tables.matrix, header=0)
        p = scatter((DATA[1, :]) , DATA[2, :],  
                xlims=(2^5, 2^12), ylims=(10^-7, 10^-2), 
                xscale=:log2, yscale=:log10,
                title="Deviation in MC phase L→∞", xlabel=L"L",
                ylabel=L"J-1/4", 
                label="α=1, β=1, p2=1", dpi=300)
        # linear fit
        @. model(x, par) = par[1]*x^-1
        xdata = DATA[1, :]
        ydata = DATA[2, :]
        par0 = [1.]
        fit = curve_fit(model, xdata, ydata, par0)
        params = @. round(fit.param, digits=2)
        println(params)
        
        # plot fit
        plot!(xdata, model(xdata, params), label=L"Fit with $ax^{-1}$")
        display("image/png", p)
end

function plot_rholeft_p2(path)
        DATA = CSV.read(path, Tables.matrix, header=0)
        p = scatter(DATA[1, :], DATA[2, :],
                xlabel=L"p_2", 
                ylabel=L"|\langle \rho_{left} \rangle - \alpha|",
                label="α=0.2, β=0.8, L=500", dpi=300)
        display("image/png", p) # export as png
end

function plot_critical_p2_fromrholeft(path)
        DATA = CSV.read(path, Tables.matrix, header=0)
        A = 0.0:0.05:0.5
        P2 = DATA[1, :]
        P2C = zeros(2, length(A))
        for (i, LEFT) in enumerate(eachrow(DATA[2:end, :]))
                # heap from the right
                P2C[1, i] = P2[(LEFT .- 0.03) .< 0][1]
                # heap from the left
                P2C[2, i] = P2[(LEFT .- 1 ./ (1 .+ LEFT) .+ 0.2 .+ 0.03) .> 0][end]
        end
        xdata = A
        ydata = vec(mean(P2C, dims=1))
        yerr = abs.((P2C[1, :] .- P2C[2, :]) ./ 2)
        p = scatter(xdata, ydata, yerr=yerr, 
                xlabel=L"\alpha", 
                ylabel=L"p_{2,c}",
                label="β=0.8, L=500", dpi=300)
        # linear fit
        @. model(x, par) = par[1]*x
        par0 = [1.0]
        fit = curve_fit(model, xdata[1:7], ydata[1:7], par0)
        params = @. round(fit.param, digits=2)

        # plot fit
        plot!(xdata[1:7], model(xdata[1:7], params), label=L"Fit with $ax$")
        plot!(xdata[7:end], model(xdata[7:end], params), ls=:dash, label=:none)
        display("image/png", p) # export as png
end

function plot_J_p2(path, α, β)
        DATA = CSV.read(path, Tables.matrix, header=0)
        println(DATA[1,1])
        p = scatter(DATA[1, :] , DATA[2, :],
                xlabel=L"p_2", 
                ylabel=L"J - \frac{1}{4}",
                label="α=$α, β=$β, L=500", dpi=300)
        vline!([α/(1-α)], lw=2, 
                label=L"p_{2,c}=\frac{\alpha}{1-\alpha}")
        x = range(0, 1, 1000)
        plot!(x, x ./ (x .+ 1).^2 .- 0.25, 
                label=L"J_{MC} - 1/4", lw=2)
        plot!(x, zeros(1000) .+ α*(1-α) .- 0.25, 
                label=L"J_{LD}-1/4", lw=2)
        display("image/png", p) # export as png
end
# plot_current_map("current-200-impurity.csv")
# plot_J_p2("J-p2.csv", 0.4, 0.8)
# plot_rholeft_p2("rholeft-p2.csv")
# plot_fs_impurity_MC_current_deviation("fs-impurity-MC-current-deviation.csv")
plot_critical_p2_fromrholeft("multiple-rholeft-p2.csv")