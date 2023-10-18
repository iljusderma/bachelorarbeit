import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_rates(string):
    return np.array(string.split('/')[-1][:-4].split('-')[:-1]).astype(float)/100

def plot_errorbars(path, average_density, density_err, average_current):
    rates = get_rates(path)
    iterations = len(average_density)
    fig = plt.figure()
    # plot curve with error band 
    ax1 = fig.add_subplot(221, title='Evolution of density', xlabel='iterations', ylabel='density')
    step = 10
    ax1.errorbar(np.arange(iterations)[::step], average_density[::step], yerr=density_err[::step], fmt=',', ecolor='grey', color='orange')

    # plot analytical prediction
    if rates[0] > rates[1] and rates[1] < 0.5:
        ax1.plot(np.full((iterations,), 1-rates[1]), color='red', label='HD')
    if rates[0] < rates[1] and rates[0] < 0.5:
        ax1.plot(np.full((iterations,), rates[0]), color='red', label='LD')
    elif rates[0] > 0.5 and rates[1] > 0.5:
        ax1.plot(np.full((iterations,), 0.5), color='red', label='MC')
    else:
        pass
    ax1.legend()
    # plot error curve
    ax2 = fig.add_subplot(222, title='Error evolution over the iterations', xlabel='iterations', ylabel='Error')
    ax2.plot(density_err)
    # plot current
    ax3 = fig.add_subplot(223, title='Evolution of the current', xlabel='iterations', ylabel='current')
    ax3.plot(average_current)
    plt.show()

# read csv
path = "/home2/ilja/Documents/Bachelorarbeit/bachelorarbeit/60-30-70-0-1000.csv"
df = pd.read_csv(path, header=None, index_col=False)
density = df.to_numpy()
# calculate the mean over all rows
average_density = np.mean(density, axis=0)
# calculate the errors / variance
# Schwankung: dx = mean(x**2) - mean(x)**2
dx = np.mean(density**2, axis=0) - average_density**2
# Standardabweichung
err = np.std(density, axis=0)

plot_errorbars(path, average_density, err)