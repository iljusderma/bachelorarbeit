import csv

with open('persons.csv', 'w', newline='') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(['alpha', 'beta', 'p', 'q', 'iterations', 'density_values'])