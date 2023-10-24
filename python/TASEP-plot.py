import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_rates(string):
    return np.array(string.split('/')[-1][:-4].split('-')[:-1]).astype(float)/100

def read_data(path):
    df = pd.read_csv(path, header=None, index_col=False)
    return df.to_numpy()

def plot_errorbars(path, average_data, err):
    rates = get_rates(path)
    iterations = len(average_data)
    fig = plt.figure()
    # plot curve with error band 
    ax1 = fig.add_subplot(111, title='Evolution over time', xlabel='iterations', ylabel='current')
    step = 1
    ax1.errorbar(np.arange(iterations)[::step], average_data[::step], yerr=err[::step], fmt='.', ms=3, ecolor='grey', color='orange')
    '''
    # plot analytical prediction
    if rates[0] > rates[1] and rates[1] < 0.5:
        ax1.plot(np.full((iterations,), 1-rates[1]), color='red', label='HD')
    if rates[0] < rates[1] and rates[0] < 0.5:
        ax1.plot(np.full((iterations,), rates[0]), color='red', label='LD')
    elif rates[0] > 0.5 and rates[1] > 0.5:
        ax1.plot(np.full((iterations,), 0.5), color='red', label='MC')
    else:
        pass'''
    ax1.legend()
    # plot error curve
    #ax2 = fig.add_subplot(122, title='Error evolution over the iterations', xlabel='iterations', ylabel='Error')
    #ax2.plot(err)
    plt.show()

# read csv
'''
current_path = "/home/ilja/Documents/coding/bachelorarbeit/60-30-70-0-1000-current.csv"
density_path = "/home/ilja/Documents/coding/bachelorarbeit/60-30-70-0-1000-density.csv"
density = read_data(density_path)
current = read_data(current_path)
# calculate the mean over all rows
average_density = np.mean(density, axis=0)
average_current = np.mean(current, axis=0)
# calculate the errors / variance
# Schwankung: dx = mean(x**2) - mean(x)**2
dx = np.mean(density**2, axis=0) - average_density**2
# Standardabweichung
density_err = np.std(density, axis=0)
current_err = np.std(current, axis=0)
# plot results
plot_errorbars(current_path, average_current, current_err)'''

path = "/home/ilja/Documents/coding/bachelorarbeit/python/current-50.csv"
CURRENT = read_data(path)
print(CURRENT.shape)
fig, ax = plt.subplots()
im = ax.imshow(CURRENT, origin='lower', cmap='plasma', extent=[0,1,0,1])
ax.set_xlabel(r"$\beta$")
ax.set_ylabel(r"$\alpha$")
fig.colorbar(im, label='Mean current', location='top')

plt.show()