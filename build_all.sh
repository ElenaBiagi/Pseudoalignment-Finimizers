folder=$1
if [ ! -d "$folder/index" ]; then
  mkdir ${folder}/index
fi

m=0.01
t=1
for i in $(seq 1 10);
do
    source /module_load.sh
    echo "0.a. Build the index for Finimizers"
    .$HOME/Pseudoalignment-Finimizers/SBWT/build/bin/sbwt build -i $HOME/${folder}/${folder}_concat.fa -o $HOME/${folder}/index/${folder}_concat.sbwt -k 31 --add-reverse-complements -d ./SBWT/temp/
    /usr/bin/time -l .$HOME/Pseudoalignment-Finimizers/benchmark build-fmin -i $HOME/${folder}/index/${folder}_31.sbwt -t 1 -o $HOME/${folder}/index/fmin_${folder}_31 -c $HOME/${folder}/${folder}_list.txt 1> $HOME/Pseudoalignment-Finimizers/results/build_${folder}-k31.stdout 2> $HOME/Pseudoalignment-Finimizers/results/build_${folder}-k31.stderr
    #--lcs $HOME/Pseudoalignment-Finimizers/results/${folder}_concat_31.LCS.sdsl 

    echo "0.b. Build the index for Themisto"
    $HOME/themisto/build/bin/themisto build -k 31 -i $HOME/${folder}/${folder}_list.txt --index-prefix $HOME/${folder}/index/themisto_${folder}_31 --temp-dir $HOME/themisto/temp --mem-gigas 2 --n-threads 4 --file-colors
    source /module_load_cobs.sh 
    echo "0.c. Build the index for COBS"
    $HOME/cobs/build/src/cobs compact-construct --clobber $HOME/${folder}/num/ ${folder}/index/cobs_${folder}.cobs_compact
    source /module_load.sh
    echo "0.d. Build the index for BLAST"
    makeblastdb -in $HOME/${folder}/seq/${folder}_num_concat.fa -parse_seqids -blastdb_version 5 -out "${folder}/index/blastdb_${folder}" -dbtype nucl
done
