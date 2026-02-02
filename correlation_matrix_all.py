
import re
import os
import json
import argparse
import sys
from collections import defaultdict
import numpy as np
from sklearn.metrics import matthews_corrcoef
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("query_n", help="Enter the query name", type=str)

def save_scatter_pdf(
    x, y, filename,
    xlabel="X",
    ylabel="Y",
    point_size=8,
    point_color="steelblue"
):
    assert len(x) == len(y)

    plt.figure(figsize=(6, 6))

    plt.scatter(x, y, s=point_size, color=point_color)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)

    plt.tight_layout()
    plt.savefig(filename, format="pdf")
    plt.close()
               
def fill_minimap_hm(query_name, max_matches):
    """
    For each genome, read its PAF file for the query_name.
    - Only consider tp:A:P lines
    - Take column 10 (1-based) as the value
    """
    minimap_res = defaultdict(int)
    tot_matches = []
    for genome in range(1, 3683):        
        paf = f"/home/biagiele/Ecoli/reads/results/{genome}/new_{genome}_{query_name}_15.paf"
        #paf = f"/home/biagiele/Ecoli/reads/results/{genome}/{genome}_{query_name}.paf"

        with open(paf) as f:
            prev_query_i = 0 
            for line in f:
                
                if "tp:A:P" not in line:     #keep only tp:A:P
                    continue
                cols = line.split()
                if len(cols) < 10:
                    continue
                #query_n = cols[0] # query_name.number
                #query_i = re.search(r"\.(\d+)$", cols[0]) # query index

                
                query_i = re.search(r'\.(\d+)$', cols[0])
                query_i = int(query_i.group(1))
                # if query_i > 10:
                #     break
                if query_i != prev_query_i:
                    prev_query_i = query_i
                    tot_matches = [0]*int(cols[1])
                start = int(cols[2])
                end = int(cols[3])
                
                # TODO read the cigar string which starts with Z
                # cigar string = cg:Z:4=10X5I2D3=...
                # Extract cigar string
                cigar = None
                # total_matches = 0
                for col in cols[12:]:
                    if col.startswith("cg:Z:"):
                        cigar = col[5:]
                        # Parse CIGAR string properly: extract all operation-length pairs
                        # For each match (=), mark positions in max_matches
                        operations = re.findall(r'(\d+)([=XID])', cigar)
                        read_pos = start  # position in the read
                        
                        for length_str, operation in operations:
                            length = int(length_str)
                            
                            if operation == '=':
                                # Mark all matching positions with bounds check
                                end_pos = min(read_pos + length, len(max_matches[int(query_i)-1]))
                                for i in range(read_pos, end_pos):
                                    max_matches[int(query_i)-1][i] = 1
                                    tot_matches[i] = 1  # mark the matches in the single query
                                read_pos += length
                                # total_matches += length
                            elif operation == 'X':
                                # Mismatch: advance position but don't mark
                                read_pos += length
                            elif operation == 'I':
                                # Insertion in query: advance query position but not reference
                                read_pos += length
                        # calculate total matches from cigar
                        # match_nums = re.findall(r'(\d+)=', cigar)
                        # total_matches = sum(int(num) for num in match_nums)
                value = sum(tot_matches) # int(cols[9])
 
                query_idx = cols[0]
                # if total_matches != value:
                #     cerr << f"Warning: total matches {total_matches} != value {value} in {paf} for read {query_idx}\n"
                minimap_res[(query_idx, genome-1)] = max(minimap_res[(query_idx, genome-1)], value)
                # minimap_res[(query_idx, genome-1)] = minimap_res[(query_idx, genome-1)] + value

                
    return minimap_res

import gzip

def fastq_gz_lengths(path):
    lengths = []
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # quality
            lengths.append(len(seq))
    return lengths

def main():
    n_cols = 3682       # genomes
    n_tools = 4         # themisto, finimap, kaminari, minimap
    t = 0             # threshold
    
    args = parser.parse_args()
    query_n = args.query_n
    print("query_name =", query_n, flush=True)

    for query_name in [query_n]:
        """ERR10496621, "ERR10496620", "ERR10496619", "ERR10496618",
            "ERR10496617", "ERR10496616", "ERR10496615", "ERR10496755",
            "ERR10496623", "ERR10496622"
        ]:
        """
        # ------------------------------------------
        # READ ALL LINES OF ALL 4 FILES FIRST
        # ------------------------------------------

        # THEMISTO --------------------------------
        print("read themisto", flush=True)
        themisto_file = f"/home/biagiele/Ecoli/reads/results/{query_name}_themisto2_res_t{t}.txt"
        with open(themisto_file) as f:
            themisto_lines = [line.strip() for line in f]

        # FINIMAP ---------------------------------
        print("read finimap", flush=True)
        finimap_file = f"/home/biagiele/Ecoli/reads/results/bases_{query_name}_10_fmin_Ecoli_{t}.txt"
        with open(finimap_file, "r") as f:
            finimap_lines = [line.strip() for line in f]

        # KAMINARI ---------------------------------        
        print("read kaminari", flush=True)
        if (t == 0):kaminari_file = f"/home/biagiele/Ecoli/reads/results/{query_name}_kaminari_res_t0.00000001_19.txt"
        else: kaminari_file = f"/home/biagiele/Ecoli/reads/results/{query_name}_kaminari_res_t{t}_19.txt"
        with open(kaminari_file, "r") as f:
            kaminari_lines = [line.strip() for line in f]

        # read query lengths
        reads_file = f"/home/biagiele/Ecoli/reads/{query_name}.fastq.gz"
        query_lens = fastq_gz_lengths(reads_file)
        
        # MINIMAP ----------------------------------
        max_matches = []
        for r in range(len(query_lens)):
            max_matches.append([0]*query_lens[r])
        print("read minimap", flush=True)
        minimap_hm = fill_minimap_hm(query_name, max_matches)
        
        
        # number of rows is min of all 4
        R = len(themisto_lines)
        # min(len(themisto_lines),
        #         len(finimap_lines),
        #         len(kaminari_lines))
        print(R)
                

        
        # ------------------------------------------
        # BUILD M FOR EACH ROW r
        # ------------------------------------------
    
        # Collect values
        
        scatter_x = []
        scatter_fx = []
        scatter_kx = []
        scatter_fy = []
        scatter_ky = []
        scatter_y = []
        
        FP_t = 0
        FP_f = 0
        FP_k = 0
        
        FN_t = 0
        FN_f = 0
        FN_k = 0
        
        TN_t = 0
        TN_f = 0
        TN_k = 0
        
        TP_t = 0
        TP_f = 0
        TP_k = 0
        
        true_m = [0] * (R * n_cols) 
        pred_t = [0] * (R * n_cols)
        pred_f = [0] * (R * n_cols)
        pred_k = [0] * (R * n_cols)
                
        M = np.zeros((n_tools, n_cols*R), dtype=float)
        L = []
        for r in range(R):
            list_minimap = []
            # print(r+1, flush = True) 

            # THEMISTO ------------------------------
            # print("Insert themisto values", flush=True)
            rec = json.loads(themisto_lines[r])

            genomes = rec["colors"]
            bases = rec["bases_covered"]
            for genome_num, value in zip(genomes, bases):
                #if genome_num < 800:
                M[0, (r*n_cols)+genome_num] = value/query_lens[r] # divide by the read length
                
            # pairs = re.findall(r'(\d+):(\d+)', themisto_lines[r])
            #for genome_num, value in pairs:
            #    genome_num, value = int(genome_num), int(value)
            #    if genome_num < n_cols:
            #        M[0, genome_num] = value """

            # FINIMAP -------------------------------
            # print("Insert finimap values", flush=True)
            pairs = re.findall(r"(\d+):(\d+)", finimap_lines[r])
            for genome_num, value in pairs:
                genome_num, value = int(genome_num), int(value)
                if genome_num < n_cols:
                    M[1, (r*n_cols) + genome_num] = value/query_lens[r] # divide by the read length

            # KAMINARI ------------------------------
            # print("Insert kaminari values", flush=True)
            pairs = re.findall(r"\((\d+),\s*(\d+)\)", kaminari_lines[r])
            for genome_num, value in pairs:
                genome_num, value = int(genome_num), int(value)
                if genome_num < n_cols:
                    M[2, (r*n_cols) + genome_num] = value/(query_lens[r]-31+1) # divide by the number of kmers
                    
        
            # MINIMAP -------------------------------
            # print("Insert minimap values", flush=True)
            for genome_num in range(n_cols):
                key = (f"{query_name}.{r+1}", genome_num)
                if key in minimap_hm:
                    value = minimap_hm[key]
                    M[3,(r*n_cols) + genome_num] = min(value/query_lens[r], 1) # divide by the read length
                    list_minimap.append([min(value/query_lens[r], 1), genome_num])
            # maximum match = denominator 
            list_minimap.sort(key=lambda x: x[0], reverse=True)
            
            minimap_threshold = sum(max_matches[r])/query_lens[r]
            
            # Remove values below threshold
            list_minimap = [res for res in list_minimap if res[0] >= minimap_threshold * t]
            
            L.append(list_minimap)

            # Remove values below the threshold
            M[3, (r * n_cols): (r * n_cols) + n_cols] [M[3, (r * n_cols): (r * n_cols) + n_cols] < (minimap_threshold * t)] = 0

            # for genome in range(n_cols):
            #     tx = M[0, (r*n_cols) + genome]   # Themisto
            #     fx = M[1, (r*n_cols) + genome]   # Finimap
            #     kx = M[2, (r*n_cols) + genome]   # Kaminari
            #     my = M[3, (r*n_cols) + genome]   # Minimap2
            
            #     if my > 0:
            #         true_m[(r * n_cols) + genome] = 1
                
            #     if tx > 0:
            #         pred_t[(r * n_cols) + genome] = 1
                
            #     if fx > 0:
            #         pred_f[(r * n_cols) + genome] = 1
                
            #     if kx > 0:
            #         pred_k[(r * n_cols) + genome] = 1
                
            #     # FALSE POSITIVES
            #     if my > 0:
            #         # False negatives
            #         if tx == 0:
            #             FN_t+=1
            #         if fx == 0:
            #             FN_f+=1
            #         if kx == 0:
            #             FN_k+=1
            #         # True positives
            #         if tx > 0:
            #             TP_t+=1
            #         if fx > 0:
            #             TP_f+=1
            #         if kx > 0:
            #             TP_k+=1
                        
            #     elif my == 0:
            #         # true negatives
            #         if tx == 0:
            #             TN_t+=1
            #         if fx == 0:
            #             TN_f+=1
            #         if kx == 0:
            #             TN_k+=1
            #         # False positives
            #         if tx > 0:
            #             FP_t+=1
            #         if fx > 0:
            #             FP_f+=1
            #         if kx > 0:
            #             FP_k+=1
                
                
            #     # scatter plot values
            #     if tx > 0 or my > 0:   # avoid empty points
            #         scatter_x.append(tx)
            #         scatter_y.append(my)
            #     if fx > 0 or my > 0:   # avoid empty points
            #         scatter_fx.append(fx)
            #         scatter_fy.append(my)
            #     if kx > 0 or my > 0:   # avoid empty points
            #         scatter_kx.append(kx)
            #         scatter_ky.append(my)
                                
    # --------------------------------------
    # CORRELATION MATRIX
    # --------------------------------------
    # print(true_m)
    # print(pred_t)
    # print(pred_f)
    # print(matthews_corrcoef(true_m, pred_t))
    # print(matthews_corrcoef(true_m, pred_f))
    # print(matthews_corrcoef(true_m, pred_k))

    # Save M[0,:] Themisto
    with open(f"/home/biagiele/Ecoli/reads/results/{query_name}_themisto2_list_t{t}.txt", "w") as f:
        i=0
        for value in M[0, :]:
            i+=1
            f.write(f"{value}\n")
            if (i > n_cols):  
                f.write(f"\n")
                i=0
    # Save M[1,:] Finimap
    with open(f"/home/biagiele/Ecoli/reads/results/{query_name}_finimap_list_t{t}.txt", "w") as f:
        i=0
        for value in M[1, :]:
            i+=1
            f.write(f"{value}\t")
            if (i > n_cols):  
                f.write(f"\n")
                i=0
    # Save M[2,:] kaminari
    with open(f"/home/biagiele/Ecoli/reads/results/{query_name}_kaminari_list_t{t}.txt", "w") as f:
        i=0
        for value in M[2, :]:
            i+=1
            f.write(f"{value}\t")
            if (i > n_cols):  
                f.write(f"\n")
                i=0
    # Save M[3,:] Minimap
    with open(f"/home/biagiele/Ecoli/reads/results/{query_name}_minimap_list_t{t}.txt", "w") as f:
        i=0
        for value in M[3, :]:
            i+=1
            f.write(f"{value}\t")
            if (i > n_cols):  
                f.write(f"\n")
                i=0
    with open(f"/home/biagiele/Ecoli/reads/results/{query_name}_minimap_norm_pairs_t{t}.txt", "w") as f:
        for read in L:
            for value in read:
                if value[0] > 0:
                    f.write(f"{value[1]}:{value[0]}\t")
            f.write(f"\n")

    # corr = np.corrcoef(M)
    
    # print(f"Correlation matrix for {query_name}")
    # print(corr)
    # print()
            
    # if scatter_x:
    #     save_scatter_pdf(
    #         scatter_x,
    #         scatter_y,
    #         "norm_{query_name}_themisto_vs_minimap2_0.8.pdf",
    #         xlabel="Themisto",
    #         ylabel="Minimap2"
    #     )
    #     print(max(scatter_x))
    #     print(max(scatter_y))
    # if scatter_fx:
    #     save_scatter_pdf(
    #         scatter_fx,
    #         scatter_fy,
    #         f"norm_{query_name}_finimap_vs_minimap2_0.8.pdf",
    #         xlabel="Finimap",
    #         ylabel="Minimap2",
    #         point_color="chocolate",#(210, 105, 30), #
    #     )
    #     print(max(scatter_fx))
    #     print(max(scatter_fy))
    # if scatter_kx:
    #     save_scatter_pdf(
    #         scatter_kx,
    #         scatter_ky,
    #         f"norm_{query_name}_kaminari_vs_minimap2_0.8.pdf",
    #         xlabel="Kaminari",
    #         ylabel="Minimap2",
    #         point_color="red",
    #     )
    #     print(max(scatter_kx))
    #     print(max(scatter_ky))
    
    
    # print()
    # print(TP_t)
    # print(TP_f)
    # print(TP_k)
    # print()
    # print(FP_t)
    # print(FP_f)
    # print(FP_k)
    # print()
    # print(FN_t)
    # print(FN_f)
    # print(FN_k)
    # print()
    # print(TN_t)
    # print(TN_f)
    # print(TN_k)
    
    # print()
    # print("False positive rate = FP/(FP+TN)")
    # print(FP_t/(FP_t+TN_t))
    # print(FP_f/(FP_f+TN_f))
    # print(FP_k/(FP_t+TN_t))
    
    # print("'%' false positives")
    # print(FP_t/(R*n_cols))
    # print(FP_f/(R*n_cols))
    # print(FP_k/(R*n_cols))
    
    # print("False negative rate = FN/(TP+FN)")
    # print(FN_t/(FN_t+TP_t))
    # print(FN_f/(FN_f+TP_f))
    # print(FN_k/(FN_t+TP_t))
    
    # print("'%' false negatives")
    # print(FN_t/(R*n_cols))
    # print(FN_f/(R*n_cols))
    # print(FN_k/(R*n_cols))
    
    # print("Accuracy")
    # print((TP_t+TN_t)/(R*n_cols))
    # print((TP_f+TN_f)/(R*n_cols))
    # print((TP_k+TN_k)/(R*n_cols))
    
    # print("Precision")
    # print(TP_t/(TP_t +FN_t))
    # print(TP_f/(TP_f +FN_f))
    # print(TP_k/(TP_k +FN_k))
    
    # print("Recall")
    # print(TP_t/(FN_t+TP_t))
    # print(TP_f/(FN_f+TP_f))
    # print(TP_k/(FN_t+TP_t))
    

if __name__ == "__main__":
    main()