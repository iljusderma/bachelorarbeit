import csv
import numpy as np
import TASEP
import pandas as pd 

probs = np.array([0.3, 0.6, 0.8, 0])
L = 1000
iterations = 10*1000
Runs = 3
density = np.zeros((Runs, iterations))

state = TASEP.init_state(L)

title = '{}-{}-{}-{}-{}.csv'.format(*(100*probs).astype(int), L)

with open(title, 'w', newline='') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #filewriter.writerow(['N', 'Density_values'])

for r in range(Runs):
    print(r)
    density[r] = TASEP.update(iterations, L, state, probs)

df = pd.DataFrame(density)
df.to_csv(title, header=False)