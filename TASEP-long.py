import csv
import numpy as np
import time
import TASEP
import pandas as pd 

probs = np.array([0.6, 0.3, 0.7, 0])
L = 1000
iterations = 50*1000
Runs = 100   
density = np.zeros((Runs, iterations))

# create CSV
title = '{}-{}-{}-{}-{}.csv'.format(*(100*probs).astype(int), L)
with open(title, 'w', newline='') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)

# perform various simulations and saving the density in a .csv file
for r in range(Runs):
    print(r)
    density[r] = TASEP.update(iterations, L, TASEP.init_state(L), probs)

df = pd.DataFrame(density)
df.to_csv(title, header=False, index=False)