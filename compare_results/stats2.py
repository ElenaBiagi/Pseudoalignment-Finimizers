import numpy as np
from sklearn.metrics import f1_score, precision_score, recall_score
from scipy.stats import wilcoxon
import sys

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
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    
    
    genomes = int(sys.argv[3]) 

    # Read all lines from file
    true = read_data(file1_path)
    pred = read_data(file2_path)

    # Calculate F1 score
    f1_themisto = []
    precision_themisto = []
    recall_themisto = []
    for t, p in zip(true, pred):
        f1 = f1_score(t, p, average="weighted")
        precision = precision_score(t, p)
        recall = recall_score(t, p)
        
        f1_themisto.append(f1)
        precision_themisto.append(precision)
        recall_themisto.append(recall)

    # print("Mean F1 score:", np.mean(f1_themisto))

    # Random classifier
    random = np.random.randint(0, 2, size=(len(true), 43))
    f1_random = [f1_score(t, r, average="weighted") for t, r in zip(true, random)]
    #print("Random F1 score:", np.mean(f1_random))

    # Wilcoxon test
    _, p_value = wilcoxon(f1_themisto, f1_random)
    #print("P-value:", p_value)
    
    print(f'{np.mean(precision_themisto):.2f}')
    print(f'{np.mean(recall_themisto):.2f}')
    print(f'{np.mean(f1_themisto):.2f}')
    print(f'{np.mean(f1_random):.2f}')
    print(p_value)
    print()
