import sys

def parse_file(file_path):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = list(map(int, line.split()))
            row_num = parts[0]
            numbers = set(parts[1:])
            data[row_num] = numbers
    return data

def calculate_precision(TP, FP):
    return TP / (TP + FP) if (TP + FP) > 0 else 0

def calculate_recall(TP, FN):
    return TP / (TP + FN) if (TP + FN) > 0 else 0

def calculate_f1_score(precision, recall):
    return 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

def calculate_accuracy(TP, TN, FP, FN):
    return (TP + TN) / (TP + TN + FP + FN) if (TP + TN + FP + FN) > 0 else 0

def calculate_specificity(TN, FP):
    return TN / (TN + FP) if (TN + FP) > 0 else 0

def compare_files(file1_data, file2_data, genomes):
    results = {}
    tot = {
        'TP': 0, 'FP': 0, 'FN': 0, 'TN': 0,
        'Precision': 0, 'Recall': 0,
        'F1 Score': 0, 'Accuracy': 0,
        'Specificity': 0
    }
    num_rows = len(file1_data)
    
    for row_num in file1_data:
        set1 = file1_data.get(row_num, set())
        set2 = file2_data.get(row_num, set())
        
        TP = len(set1 & set2)  # True Positives
        FP = len(set2 - set1)  # False Positives
        FN = len(set1 - set2)  # False Negatives
        TN = genomes - len(set1 | set2)  # True Negatives
        
        precision = calculate_precision(TP, FP)
        recall = calculate_recall(TP, FN)
        f1_score = calculate_f1_score(precision, recall)
        accuracy = calculate_accuracy(TP, TN, FP, FN)
        specificity = calculate_specificity(TN, FP)
        
        results[row_num] = {
            'TP': TP, 'FP': FP, 'FN': FN, 'TN': TN,
            'Precision': precision, 'Recall': recall,
            'F1 Score': f1_score, 'Accuracy': accuracy,
            'Specificity': specificity
        }
        tot['TP'] += TP
        tot['FP'] += FP
        tot['FN'] += FN
        tot['TN'] += TN
        tot['Precision'] += precision
        tot['Recall'] += recall
        tot['F1 Score'] += f1_score
        tot['Accuracy'] += accuracy
        tot['Specificity'] += specificity
    
    avg_tot = {key: value / num_rows for key, value in tot.items()}
    
    return results, avg_tot

def main():
    file1_path = sys.argv[1] # THEMISTO
    file2_path = sys.argv[2] # Finimizers
    
    genomes = int(sys.argv[3]) # Numbers range from 0 to 42 inclusive

    file1_data = parse_file(file1_path)
    file2_data = parse_file(file2_path)

    comparison_results, avg_tot = compare_files(file1_data, file2_data, genomes)

    print(f'{avg_tot['TP']:.2f}')
    print(f'{avg_tot['FP']:.2f}')
    print(f'{avg_tot['FN']:.2f}')
    print(f'{avg_tot['TN']:.2f}')
    print(f'{avg_tot['Accuracy']:.2f}')
    print(f'{avg_tot['Specificity']:.2f}')
    #print(f'{avg_tot['Precision']:.2f}')
    #print(f'{avg_tot['Recall']:.2f}')
    #print(f'{avg_tot['F1 Score']:.2f}')
    

if __name__ == "__main__":
    main()