import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_rates(string):
    return np.array(string.split('/')[-1][:-4].split('-')[:-1]).astype(float)/100

# read csv
path = "/home2/ilja/Documents/Bachelorarbeit/bachelorarbeit/60-30-70-0-1000.csv"
df = pd.read_csv(path, header=None, index_col=False)
density = df.to_numpy()
# calculate the mean over all rows
average_density = np.mean(density, axis=0)
rates = get_rates(path)
print(rates)
iterations = 50*1000

# plot
fig = plt.figure()
ax = fig.add_subplot(111, title='Evolution of density', xlabel='iterations', ylabel='density')
ax.plot(np.arange(len(average_density)), average_density, marker=',', lw=0)
if rates[0] > rates[1] and rates[1] < 0.5:
    ax.plot(np.full((iterations,), 1-rates[1]), color='red', label='HD')
if rates[0] < rates[1] and rates[0] < 0.5:
    ax.plot(np.full((iterations,), rates[0]), color='red', label='LD')
elif rates[0] > 0.5 and rates[1] > 0.5:
    ax.plot(np.full((iterations,), 0.5), color='red', label='MC')
else:
    pass
ax.legend()
plt.show()