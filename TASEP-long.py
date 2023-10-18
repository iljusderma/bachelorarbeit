import csv
import numpy as np
import TASEP
import pandas as pd

def save_data(data, name):
    # create CSV
    title = '{}-{}-{}-{}-{}-' + name + '.csv'.format(*(100*rates).astype(int), L)
    with open(title, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
    # save data to CSV
    data_df = pd.DataFrame(data)
    data_df.to_csv(title, header=False, index=False)

rates = np.array([0.6, 0.3, 0.7, 0])
L = 1000
iterations = 50*1000
Runs = 100
current_stepsize = 100
density = np.zeros((Runs, iterations))
current = np.zeros((Runs, (iterations - 1)//current_stepsize))

# perform various simulations
for r in range(Runs):
    print(str(round(r/Runs, 2)) + '%')
    chain = TASEP.Chain(L, rates)
    chain.initialize_state()
    density[r], current[r] = TASEP.iterate(iterations, chain, current_stepsize)

save_data(density, [ i for i, j in locals().items() if j == density][0])
save_data(current, [ i for i, j in locals().items() if j == current][0])