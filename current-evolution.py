import numpy as np
import TASEP
import matplotlib.pyplot as plt
import pandas as pd
import csv

def save_data(data, name):
    # create CSV
    title = name + '.csv'
    with open(title, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
    # save data to CSV
    data_df = pd.DataFrame(data)
    data_df.to_csv(title, header=False, index=False)

rates = np.array([0.6, 0.3, 1, 0])
ALPHA, BETA = np.linspace(0, 1, 10), np.linspace(0, 1, 10)
CURRENT = np.zeros((10, 10))
L = 500
iterations = 10*1000
current_stepsize = 50

for index, current in np.ndenumerate(CURRENT):
    print(index)
    a, b = index
    rates[0], rates[1] = ALPHA[a], BETA[b]
    chain = TASEP.Chain(L, rates)
    chain.initialize_state()
    CURRENT[index] = np.mean(TASEP.iterate(iterations, chain, current_stepsize)[-1])

fig = plt.figure()
ax = fig.add_subplot(111, aspect=1.0, title='', xlabel='', ylabel='')
ax.imshow(CURRENT)

save_data(CURRENT, 'current')