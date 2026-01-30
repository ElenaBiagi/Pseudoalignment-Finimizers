
import re
import os
import json
import argparse
import sys
from collections import defaultdict
import numpy as np
import gzip
from sklearn.metrics import matthews_corrcoef, roc_curve, auc, precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("query_n", help="Enter the query name", type=str)
             
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
    tt = 0         # threshold
    print("Themisto Threshold:", tt)
    args = parser.parse_args()
    query_n = args.query_n
    print("query_name =", query_n, flush=True)

    for query_name in [query_n]:
        """ERR10496621, "ERR10496620", "ERR10496619", "ERR10496618",
            "ERR10496617", "ERR10496616", "ERR10496615", "ERR10496755",
            "ERR10496623", "ERR10496622"
        """
        # ------------------------------------------
        # READ ALL LINES OF ALL 4 FILES FIRST
        # ------------------------------------------

        # THEMISTO --------------------------------
        # print("read themisto", flush=True)
        themisto_file = f"/home/biagiele/Ecoli/reads/results/{query_name}_themisto2_res_t{tt}.txt"
        with open(themisto_file) as f:
            themisto_lines = [line.strip() for line in f]

        tf = 0
        print("F & K Threshold:", tf)
        # FINIMAP ---------------------------------
        # print("read finimap", flush=True)
        finimap_file = f"/home/biagiele/Ecoli/reads/results/bases_{query_name}_10_fmin_Ecoli_{tf}.txt"
        with open(finimap_file, "r") as f:
            finimap_lines = [line.strip() for line in f]

        # KAMINARI ---------------------------------        
        # print("read kaminari", flush=True)
        if (tf == 0):
            kaminari_file = f"/home/biagiele/Ecoli/reads/results/{query_name}_kaminari_res_t0.00000001_19.txt"
        else:
            kaminari_file = f"/home/biagiele/Ecoli/reads/results/{query_name}_kaminari_res_t{tf}_19.txt"
        with open(kaminari_file, "r") as f:
            kaminari_lines = [line.strip() for line in f]

        # read query lengths
        reads_file = f"/home/biagiele/Ecoli/reads/{query_name}.fastq.gz"
        query_lens = fastq_gz_lengths(reads_file)
        
        # MINIMAP ----------------------------------
        max_matches = []
        for r in range(len(query_lens)):
            max_matches.append([0]*query_lens[r])
        # print("read minimap", flush=True)
        minimap_hm = fill_minimap_hm(query_name, max_matches)
        
        
        # number of rows is min of all 4
        R = len(themisto_lines)
        # min(len(themisto_lines),
        #         len(finimap_lines),
        #         len(kaminari_lines))
        print(R*n_cols)
        
        # ------------------------------------------
        # BUILD M FOR EACH ROW r
        # ------------------------------------------
    
        # Collect values
        true_m = [0] * (R * n_cols) 
        pred_t = [0] * (R * n_cols)
        pred_f = [0] * (R * n_cols)
        pred_k = [0] * (R * n_cols)
                
        M = np.zeros((n_tools, n_cols*R), dtype=float)

        for r in range(R):
            # print(r+1, flush = True) 

            # THEMISTO ------------------------------
            # print("Insert themisto values", flush=True)
            rec = json.loads(themisto_lines[r])

            genomes = rec["colors"]
            bases = rec["bases_covered"]
            for genome_num, value in zip(genomes, bases):
                #if genome_num < 800:
                M[0, (r*n_cols)+genome_num] = value/query_lens[r] # divide by the read length

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
            tm = 0
            if r == 0: 
                print("Minimap threshold:", tm)
            for genome_num in range(n_cols):
                key = (f"{query_name}.{r+1}", genome_num)
                if key in minimap_hm:
                    value = minimap_hm[key]
                    M[3,(r*n_cols) + genome_num] = min(value/query_lens[r],1) # divide by the read length

            minimap_threshold = sum(max_matches[r])/query_lens[r]

            # Remove values below the threshold
            # M[3, (r * n_cols): (r * n_cols) + n_cols] [M[3, (r * n_cols): (r * n_cols) + n_cols] < (minimap_threshold * t)] = 0
            # # Set everything else to 1
            # M[3, (r * n_cols): (r * n_cols) + n_cols] [M[3, (r * n_cols): (r * n_cols) + n_cols] != 0] = 1
            start = r * n_cols
            end = start + n_cols

            if (tm>0): M[3, start:end] = (M[3, start:end] >= (minimap_threshold * tm)).astype(int)
            else: M[3, start:end] = (M[3, start:end] > 0).astype(int)
    
    pred_t = (M[0, :] > 0).astype(int)
    pred_f = (M[1, :] > 0).astype(int)
    pred_k = (M[2, :] > 0).astype(int)
    
    labels = M[3, :] # Minimap2 
    
    tools = {
        "Themisto": M[0, :],
        "Finimap ": M[1, :],
        "Kaminari": M[2, :],
    }

    plt.figure(figsize=(6, 6))
    
    print(f"#zero:{(labels==0).sum()}, #ones: {(labels==1).sum()}")

    for name, scores in tools.items():
        fpr, tpr, _ = roc_curve(labels, scores)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"{name} (AUC={roc_auc:.3f})")

    plt.plot([0, 1], [0, 1], "k--", alpha=0.6)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curves")
    plt.legend(loc="lower right")
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(f"roc_all_tools_{query_n}_{tt}-{tf}-{tm}.pdf")
    plt.close()
    
    
    plt.figure(figsize=(6, 6))

    for name, scores in tools.items():
        precision, recall, _ = precision_recall_curve(labels, scores)
        ap = average_precision_score(labels, scores)
        plt.plot(recall, precision, label=f"{name} (AP={ap:.3f})")

    # Baseline = positive prevalence
    baseline = labels.mean()
    plt.hlines(baseline, 0, 1, linestyles="dashed", colors="gray",
            label=f"Baseline ({baseline:.3f})")

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision–Recall Curves")
    plt.legend(loc="upper right")
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(f"pr_all_tools_{query_n}_{tt}-{tf}-{tm}.pdf")
    plt.close()
    
    print("Positive rate:", labels.mean())
    print("Minimap:", labels.mean() * R * n_cols)
    for name, scores in tools.items():
        print(name, np.count_nonzero(scores))
    
    # --------------------------------------
    # CORRELATION MATRIX
    # --------------------------------------
    print(matthews_corrcoef(labels, pred_t))
    print(matthews_corrcoef(labels, pred_f))
    print(matthews_corrcoef(labels, pred_k))

    corr = np.corrcoef(M)
    
    print(f"Correlation matrix for {query_name}")
    print(corr)
    print()
    
if __name__ == "__main__":
    main()