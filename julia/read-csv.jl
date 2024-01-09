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
        h = heatmap(ALPHA, BETA, FLUX, 
                topmargin=7mm, 
                xlabel=L"Exit rate $\beta$", 
                ylabel=L"Entry rate $\alpha$", 
                aspectratio=true, xlims=(0,1), 
                legendfont=14)
        contour!(ALPHA, BETA, Z, levels=range(0, p2/(1+p2)^2, 6)[2:end], color=:lightrainbow, lw=3)
        hline!([p2/(1+p2)], ls=:dash, color=:black, label=L"\alpha_c", legend=:topright)
        vline!([p2/(1+p2)], ls=:dash, color=:black, label=L"\beta_c")
        annotate!(1.1, 1.05, (L"Current $J$", 11))
        savefig("plot.pdf")
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
                title="α=0.8, β=0.8, d=0.3", 
                titleloc=:left, 
                ylims=(0, 0.06), xlabel=L"L",
                ylabel=L"\rho_{\mathrm{left, approx}} - \langle \rho_{\mathrm{left}} \rangle", 
                label=false, 
                legendfont=12)
        deviation = mean(DATA[2, :])
        deviation_err = std(DATA[2, :])/sqrt(length(DATA[2, :]))
        println(deviation, deviation_err)
        hline!([deviation], yerr=[deviation_err], label="average", lw=2)
        savefig("plot.pdf")
end

function plot_fs_impurity_MC_current_deviation(path)
        DATA = CSV.read(path, Tables.matrix, header=0)
        p = scatter((DATA[1, :]) , DATA[2, :], yerr=DATA[3, :], 
                xlims=(2^4, 2^12), ylims=(10^-7, 10^-2), 
                xscale=:log2, yscale=:log10,
                title="α=1, β=1, d=1", 
                titleloc=:left, 
                xlabel=L"L", ylabel=L"J-1/4", 
                label=false, rightmargin=3mm)
        # linear fit
        @. model(x, par) = par[1]*x^-1
        xdata = DATA[1, :]
        ydata = DATA[2, :]
        par0 = [1.]
        fit = curve_fit(model, xdata, ydata, par0)
        params = @. round(fit.param, digits=2)
        println(params)
        
        # plot fit
        plot!(xdata, model(xdata, params), label=L"Fit with $aL^{-1}$", legendfont=12)
        # xticks
        xs = 2 .^ (5:1:9)
        xticks!(p, xs, string.(xs)) 
        savefig("plot.pdf")
end

function plot_V_p(path, transition)
        p, RHO10, RHO500 = CSV.read(path, Tables.matrix, header=0)
        gridsize = length(p)
        scatter(p, RHO10, label="L=10", 
                xlabel=L"p", ylabel=L"V", 
                msw=0, ylims=[0, 1])
        scatter!(p, RHO500, label="L=500", msw=0)
        # draw limit L → ∞
        y = zeros(gridsize)
        if transition == 1
                α0 = 0.1
                mid = Int(floor(0.5*length(y)))
                y[1:mid] .= -0.5 .*p[1:mid] .+ α0
                y[mid+1:end] .= -0.5 .*p[mid+1:end] .+ (1-α0)
        else
                mid = Int(floor(0.5*length(y)))
                y[1:mid] .= 1 .- (p[1:mid] .+ ALPHA[1:mid]) # ∼1-β
                y[mid+1:end] .= 0.5
        end
        plot!(p, y, label="L ⟶ ∞")
        savefig("plot.pdf")
end

function plot_rholeft_d(path, α)
        DATA = CSV.read(path, Tables.matrix, header=0)
        p = scatter(DATA[1, :], (DATA[2, :]),
                title="α=$α, β=0.8, L=200",
                titleloc=:left, 
                xlabel=L"d", 
                ylabel=L"|\langle \rho_{\mathrm{left}} \rangle - \alpha|",
                label=false,
                legendfont=12)
        d = DATA[1, :]
        index=length(d[d .< 0.23])
        plot!(DATA[1, 1:index], abs.(1 ./(DATA[1, 1:index] .+ 1) .- α),
        lw=2, label=false)
        plot!(DATA[1, index:index+1], [1 ./(DATA[1,index] .+1) .- α, 0],
        lw=2, ls=:dash, label=false, color=palette(:default)[2])
        plot!(DATA[1, index+1:end], zeros(length(DATA[1, index+1:end])),
        lw=2, label=false, color=palette(:default)[2])
        annotate!(0.25, 0.75, (L"|\frac{1}{1+d} - \alpha|", 12))
        annotate!(0.65, 0.05, (L"|\alpha - \alpha|", 12))
        savefig("plot.pdf")
end

function plot_critical_d_fromrholeft(path)
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
                title="β=0.8, L=500", 
                titleloc=:left, 
                xlabel=L"\alpha", 
                ylabel=L"d_c",
                label=false, legendfont=12, 
                legend=:outerright)
        # linear fit
        @. model(x, par) = par[1]*x
        par0 = [1.0]
        fit = curve_fit(model, xdata[1:7], ydata[1:7], par0)
        params = @. round(fit.param, digits=2)

        # plot fit
        plot!(xdata[1:7], model(xdata[1:7], params), label=L"Fit with $ax$", lw=2)
        plot!(xdata[7:end], model(xdata[7:end], params), ls=:dash, label=:none, lw=2)
        α = 0:0.005:0.5
        plot!(α, α ./ (1 .- α), label=L"\frac{\alpha}{1-\alpha}", lw=2)
        savefig("plot.pdf")
end

function plot_J_d(path, α, β)
        DATA = CSV.read(path, Tables.matrix, header=0)
        println(DATA[1,1])
        p = scatter(DATA[1, :] , DATA[2, :],
                title="α=$α, β=$β, L=500", titleloc=:left, 
                xlabel=L"d", 
                ylabel=L"J - \frac{1}{4}",
                label=false, 
                legend=:bottomright, legendfont=12)
        vline!([α/(1-α)], lw=2, 
                label=L"$d_c$ = $\frac{\alpha}{1-\alpha}$")
        x = range(0, 1, 1000)
        plot!(x, x ./ (x .+ 1).^2 .- 0.25, 
                label=L"J_{MC} - 1/4", lw=2)
        plot!(x, zeros(1000) .+ α*(1-α) .- 0.25, 
                label=L"J_{LD}-1/4", lw=2)
        savefig("plot.pdf")
end

# plot_current_map("current-200-impurity.csv")
# plot_V_p("V-p-1order.csv", 1)
# plot_J_d("J-d-0208-500.csv", 0.2, 0.8)
# plot_rholeft_d("julia/critical-from-rholeft/rholeft-d-02.csv", 0.2)
# plot_fs_impurity_MC_current_deviation("fs-impurity-MC-current-deviation.csv")
plot_fs_impurity_MC_density_deviation("fs-impurity-MC-density-deviation.csv")