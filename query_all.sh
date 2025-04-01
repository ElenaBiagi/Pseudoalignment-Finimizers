folder=$1
if [ ! -d "$folder/results" ]; then
  mkdir ${folder}/results
fi

m=0.01
t=1
for i in $(seq 1 10);
do
    source /module_load.sh

    echo "3.a. Pseudoalignment with Finimizers"

    /usr/bin/time -l $HOME/Pseudoalignment-Finimizers/benchmark search-fmin -i $HOME/${folder}/index/fmin_${folder}_31 -q $HOME/${folder}/queries/queries_${folder}_1000_200_${i}.fa -t $t -o $HOME/${folder}/results/out_fi_1000_200_${i}_${m}.txt 1> $HOME/Pseudoalignment-Finimizers/results/search_${folder}-unitigs-k31.stdout 2> $HOME/Pseudoalignment-Finimizers/results/search_${folder}-unitigs-k31.stderr
    /usr/bin/time -l $HOME/Pseudoalignment-Finimizers/benchmark search-fmin -i $HOME/${folder}/index/fmin_${folder}_31 -q $HOME/${folder}/queries/m_queries_${folder}_1000_200_${i}_${m}.fa -t $t -o $HOME/${folder}/results/m_out_fi_1000_200_${i}_${m}.txt 1> $HOME/Pseudoalignment-Finimizers/results/m_search_${folder}-unitigs-k31.stdout 2> $HOME/Pseudoalignment-Finimizers/results/m_search_${folder}-unitigs-k31.stderr


    echo "3.b. Pseudoalignment with Themisto"

    $HOME/themisto/target/release/themisto2 dump-pseudoalignment-data -i $HOME/${folder}/index/themisto_${folder}_31.thm -q $HOME/${folder}/queries/m_queries_${folder}_1000_200_${i}_${m}.fa -m 0 > $HOME/${folder}/results/m_out_t_1000_200_${i}_${m}.txt
    python3 $HOME/parse_json_thm2.py $HOME/${folder}/results/m_out_t_1000_200_${i}_${m}.txt > $HOME/${folder}/results/m_out_ti_1000_200_${i}_${m}.txt

    $HOME/themisto/target/release/themisto2 dump-pseudoalignment-data -i $HOME/${folder}/index/themisto_${folder}_31.thm -q $HOME/${folder}/queries/queries_${folder}_1000_200_${i}.fa -m 0 > $HOME/${folder}/results/out_t_1000_200_${i}.txt
    python3 $HOME/parse_json_thm2.py $HOME/${folder}/results/out_t_1000_200_${i}.txt > $HOME/${folder}/results/out_ti_1000_200_${i}.txt


    echo "3. Local alignment with BLAST"
    blastn -db $HOME/${folder}/index/blastdb_${folder} -query $HOME/${folder}/queries/queries_${folder}_1000_200_${i}.fa -out $HOME/${folder}/results/out_b_1000_200_${i}_$m.txt
    python3 $HOME/parse_blast.py $HOME/${folder}/results/out_b_1000_200_${i}_${m}.txt $HOME/${folder}/results/out_bi_1000_200_${i}_${m}.txt

    blastn -db $HOME/${folder}/index/blastdb_${folder} -query $HOME/${folder}/queries/m_queries_${folder}_1000_200_${i}_${m}.fa -out $HOME/${folder}/results/m_out_b_1000_200_${i}_$m.txt
    python3 $HOME/parse_blast.py $HOME/${folder}/results/m_out_b_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_bi_1000_200_${i}_${m}.txt

    source /module_load_cobs.sh
    echo "3.c. Pseudoalignment with COBS"
    $HOME/cobs/build/src/cobs query -i $HOME/${folder}/index/cobs_${folder} -f $HOME/${folder}/queries/queries_${folder}_1000_200_${i}.fa > $HOME/${folder}/results/c_out${i}.txt
    python3 $HOME/cobs/parse_cobs_output.py $HOME/${folder}/results/c_out${i}.txt $HOME/${folder}/results/ci_out${i}.txt

    $HOME/cobs/build/src/cobs query -i $HOME/${folder}/index/cobs_${folder} -f $HOME/${folder}/queries/m_queries_${folder}_1000_200_${i}_${m}.fa > $HOME/${folder}/results/m_c_out${i}.txt
    python3 $HOME/cobs/parse_cobs_output.py $HOME/${folder}/results/m_c_out${i}.txt $HOME/${folder}/results/m_ci_out${i}.txt

done


