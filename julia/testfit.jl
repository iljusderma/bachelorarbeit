@. model(x, p) = p[1]*exp(p[2]*(x-p[3])) + p[4]
t = 100.0
xdata = 1.0:t
ydata = range(0, 1, length=trunc(Int, t))
p0 = [0.5, 0.5, 0.5, 0.5]

fit = curve_fit(model, xdata, ydata, p0)
params = fit.param