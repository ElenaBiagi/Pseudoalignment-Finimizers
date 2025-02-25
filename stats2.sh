import sys
import numpy as np
from sklearn.metrics import precision, recall, f1_score
from scipy.stats import wilcoxon

def read_data(file_name):
    # Read all lines from file
    data = []
    with open(file_name, "r") as file:
        for line in file:
            d = np.zeros((43))
            for i in line.strip().split()[1:]:
                d[int(i)] = 1
            data.append(d)
    return np.array(data)

if __name__ == "__main__":
    # Read all lines from file
    true = sys.argv[1] # THEMISTO
    pred = sys.argv[2] # Finimizers
    
    genomes = int(sys.argv[3]) # Numbers range from 0 to 42 inclusive


    # Calculate F1 score, precision and recall
    f1_themisto = []
    precision_themisto = []
    recall_themisto = []
    for t, p in zip(true, pred):
        f1 = f1_score(t, p, average="weighted")
        f1_themisto.append(f1)
        precision = precision_score(t, p)
        precision_themisto.append(precision)
        recall = recall_score(t, p)
        recall_themisto.append(recall)

    #print(np.mean(f1_themisto))

    # Random classifier
    random = np.random.randint(0, 2, size=(len(true), genomes))
    f1_random = [f1_score(t, r, average="weighted") for t, r in zip(true, random)]
    #print( np.mean(f1_random))

    # Wilcoxon test
    _, p_value = wilcoxon(f1_themisto, f1_random)
    
    print(f'{np.mean(precision_themisto):.2f}')
    print(f'{np.mean(recall_themisto):.2f}')
    print(f'{np.mean(f1_themisto):.2f}')
    print(f'{np.mean(f1_random):.2f}')
    print(p_value)
    print()