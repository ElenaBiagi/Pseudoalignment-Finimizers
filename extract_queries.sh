#!/bin/bash
#SBATCH --job-name=queries_ATB.job
#SBATCH -p short
#SBATCH -M kale
#SBATCH --chdir=/home/biagiele/Pseudoalignment-Finimizers/
#SBATCH --mem=370G
#SBATCH -t05:00:00
#SBATCH -n1
#SBATCH -c1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=elena.biagi@helsinki.fi


name=$1
folder=$2
m=0.01

if [ ! -d "$HOME/$folder/queries" ]; then
    mkdir $HOME/${folder}/queries
fi  

# for i in $(seq 6 10);
# do
#     # echo "1. Sample positive queries ($i)"
#     # #$HOME/Pseudoalignment-Finimizers/sample_random $HOME/${folder}/${name}_concat.fa 10000 1000 $i > $HOME/${folder}/queries/queries_${name}_10000_1000_${i}.fa 

#     #$HOME/Pseudoalignment-Finimizers/sample_random $HOME/${folder}/seq/${name}_num_concat.fna 10000 1000 $i > $HOME/${folder}/queries/new_queries_${name}_10000_1000_${i}.fa 
#     #$HOME/Pseudoalignment-Finimizers/sample_random $HOME/norm-hybrid-files/assemblies/hybrid/seq/TE_num_concat_blast.fna 10000 1000 $i > $HOME/${folder}/queries/queries_TE_10000_1000_${i}.fa 

#     #echo "2. Mutate queries"
#     #$HOME/Pseudoalignment-Finimizers/mutate_seq $HOME/queries/queries_${folder}_1000_200_${i}.fa $HOME/${folder}/queries/m_queries_${folder}_1000_200_${i}_ $m $i
# done

python3 /home/biagiele/Pseudoalignment-Finimizers/extract_pos_queries.py 


