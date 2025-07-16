import sys

#import rbo

from sigir2024_rbo import *

import numpy as np

def parse_file_t(file_path):
    data = {}
    
    r=0
    with open(file_path, 'r') as file:
        for line in file:
            before_colon = []
            after_colon = []
            
            parts = line.strip().split()  # Split line by spaces
            row_num = r #int(parts[0])  # First number is the row index
            r += 1
            
            for pair in parts[1:]:  # Process the remaining elements
                before, after = map(int, pair.split(":"))
                before_colon.append(str(before))
                after_colon.append(after)
            
            data[row_num] = (before_colon, after_colon)  # Store as tuple of lists
    
    return data

def parse_file(file_path):
    data = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            before_colon = []
            after_colon = []
            
            parts = line.strip().split()  # Split line by spaces
            row_num = int(parts[0])  # First number is the row index
            
            for pair in parts[1:]:  # Process the remaining elements
                before, after = map(int, pair.split(":"))
                before_colon.append(str(before))
                after_colon.append(after)
            
            data[row_num] = (before_colon, after_colon)  # Store as tuple of lists
    
    return data

def parse_file_2(file_path):
    data = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            before_colon = []
            after_colon = []
            
            parts = line.strip().split()  # Split line by spaces
            row_num = int(parts[0])  # First number is the row index
            
            for pair in parts[1:]:  # Process the remaining elements
                before, after = map(int, pair.split(":"))
                before_colon.append(before)
                after_colon.append(str(after))
            
            data[row_num] = (before_colon, after_colon)  # Store as tuple of lists
    
    return data

def parse_file_num(file_path):
    data = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            before_colon = []
            after_colon = []
            
            parts = line.strip().split()  # Split line by spaces
            row_num = int(parts[0])  # First number is the row index
            
            for pair in parts[1:]:  # Process the remaining elements
                before, after = map(int, pair.split(":"))
                before_colon.append(before)
                #after_colon.append(after)
            
            data[row_num] = before_colon
    
    return data

"""
S = [1, 2, 3]
T = [1, 3, 2]

#rbo.RankingSimilarity(S, T).rbo()
 S=from_string('red (blue green) yellow pink')
print(S)
 """


# compares genomes ids as strings 
def compare_files(file1_path, file2_path, genomes):
    file1_data = parse_file_t(file1_path) #themisto
    file2_data = parse_file(file2_path)
        
    num_rows = len(file1_data)
    
    tot_rbo=[]
    for row_num in file1_data:
        color_id1, k_matches1 = file1_data[row_num]
        color_id2, k_matches2 = file2_data[row_num]

        if color_id1==[]:
            if color_id2==[]:
                tot_rbo.append(1)
            else:
                tot_rbo.append(0)
        else:        
            rlist1 = extract_ranking(color_id1, k_matches1)
            rlist2 = extract_ranking(color_id2, k_matches2)
            
            tot_rbo.append(rbo_s(rlist1, rlist2, p=0.95, ties = 'w', score = ('ext'))) # rbo_s strings
    
    results = np.mean(tot_rbo)
    return results

def compare_matches_2(file1_path, file2_path, genomes):
    file1_data = parse_file_2(file1_path)
    file2_data = parse_file_2(file2_path)
        
    num_rows = len(file1_data)
    
    tot_rbo=[]
    for row_num in file1_data:
        color_id1, k_matches1 = file1_data[row_num]
        color_id2, k_matches2 = file2_data[row_num]

        #rlist1 = from_string(k_matches1)
        #rlist2 = from_string(k_matches2)
        rlist1 = extract_ranking(["170","170","148"], [2,3,0])
        rlist2 = extract_ranking(["170","170","148"], [2,3,0])

        #rlist1 = extract_ranking(k_matches1, color_id1)

        #rlist2 = extract_ranking(k_matches2, color_id2)
        
        
        tot_rbo.append(rbo_s(rlist1, rlist2, p=0.95, ties = 'w', score = ('ext'))) # rbo_s strings
    
    results = np.mean(tot_rbo)
    return results

# compares genome ids as numbers
def compare_genomeIDs(file1_path, file2_path, genomes):
    file1_data = parse_file_num(file1_path)
    file2_data = parse_file_num(file2_path)
        
    num_rows = len(file1_data)
    
    tot_rbo=[]
    for row_num in file1_data:
        color_id1 = file1_data[row_num]
        color_id2 = file2_data[row_num]
        
        tot_rbo.append(rbo.RankingSimilarity(color_id1, color_id2).rbo()) 
    
    results = np.mean(tot_rbo)
    return results


def main():
    # test 
    """ S = [1, 2, 3]
    T = [1, 3, 2]

    print(rbo.RankingSimilarity(S, T).rbo())
    
    S= extract_ranking(["17", "18", '25', "33", "34", "39", "0", "1", "2"],
                    [170, 170, 170, 170, 170, 170, 139, 139, 139])
    print(S) """
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    #file1_path = 'out_240.txt'  
    #file2_path = 'Salmonella_unitigs_31.stats_t240'
    
    genomes = int(sys.argv[3])
    
    results = compare_files(file1_path, file2_path, genomes) # compares genomeIDs with ties
    
    #results = compare_genomeIDs(file1_path, file2_path, genomes) # compares genomeIDs without ties
    
    #results = compare_matches_2(file1_path, file2_path, genomes) # compares #matches without ties # DOES NOT WORK

    print(results)
if __name__ == "__main__":
    main()
