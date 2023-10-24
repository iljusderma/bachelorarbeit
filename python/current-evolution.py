import numpy as np
import TASEP
import matplotlib.pyplot as plt
import pandas as pd
import csv
from pathlib import Path

def save_data(data, name):
    # create CSV
    title = name + '.csv'
    with open(title, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
    # save data to CSV
    data_df = pd.DataFrame(data)
    data_df.to_csv(title, header=False, index=False)

def read_data(path):
    df = pd.read_csv(path, header=None, index_col=False)
    return df.to_numpy()

grid_size = 50
rates = np.array([0.6, 0.3, 1, 0])
ALPHA, BETA = np.linspace(0, 1, grid_size), np.linspace(0, 1, grid_size)
CURRENT = np.zeros((grid_size, grid_size))
L = 500
iterations = 10*1000
current_stepsize = 50

# check for existing data
name = 'current-' + str(grid_size) + '.csv'
if Path(name).is_file() == True:
    existing_data = read_data(name)
    print('read')
    if len(existing_data[0]) == grid_size:
        CURRENT[0:len(existing_data[:,0])] = existing_data

for index, current in np.ndenumerate(CURRENT):
    if current != 0:
         pass
    else:
        print(index)
        a, b = index
        rates[0], rates[1] = ALPHA[a], BETA[b]
        chain = TASEP.Chain(L, rates)
        chain.initialize_state()
        CURRENT[index] = np.mean(TASEP.iterate(iterations, chain, current_stepsize)[-1])
        save_data(CURRENT, 'current-' + str(grid_size))

fig, ax = plt.subplots()
im = ax.imshow(CURRENT)
plt.show()