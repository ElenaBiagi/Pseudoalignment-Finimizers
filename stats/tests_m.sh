## Build the index 
# Salmonella   
echo "0.a Build the index for Finimizers"

#ok
#./SBWT/build/bin/sbwt build -i ../Salmonella/Salmonella_concat.fa -o ../Salmonella/Salmonella_concat.sbwt -k 31 --add-reverse-complements -d ./SBWT/temp/
#/usr/bin/time -l ./benchmark build-fmin -u /Users/biagiele/try/Salmonella/Salmonella_concat.fa -i /Users/biagiele/try/Salmonella/Salmonella_concat.sbwt --lcs /Users/biagiele/try/Pseudoalignment-Finimizers/results/Salmonella_concat_31.LCS.sdsl -t 1 -o /Users/biagiele/try/Pseudoalignment-Finimizers/results/Salmonella_concat_31 -c /Users/biagiele/try/Salmonella/Salmonella_list.txt 1> /Users/biagiele/try/Pseudoalignment-Finimizers/results/build_Salmonella-k31.stdout 2> /Users/biagiele/try/Pseudoalignment-Finimizers/results/build_Salmonella-k31.stderr

echo "0.b Build the index for Themisto"
#/Users/biagiele/try/themisto/build/bin/themisto build -k 31 -i /Users/biagiele/try/Pseudoalignment-Finimizers/Salmonella/Salmonella_list.txt --index-prefix /Users/biagiele/try/themisto/my_index --temp-dir /Users/biagiele/try/themisto/temp --mem-gigas 2 --n-threads 4 --file-colors

rm /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats_1000_200.txt

m=0.01
t=1
for i in $(seq 1 10);
do
    # sample 100 strings of length 200
    echo "1. Sample positive queries"
    #/Users/biagiele/try/Pseudoalignment-Finimizers/sample_random /Users/biagiele/try/Pseudoalignment-Finimizers/Salmonella_concat.fa 100 200 > /Users/biagiele/try/Pseudoalignment-Finimizers/Salmonella/queries/queries_Salmonella_1000_200_${i}.fa 
    ##/Users/biagiele/try/Pseudoalignment-Finimizers/sample_random /Users/biagiele/try/Salmonella/Salmonella_concat.fa 1000 200 > /Users/biagiele/try/queries/queries_Salmonella_1000_200_${i}.fa 

    echo "1.2. Mutate queries"

    ##/Users/biagiele/try/Pseudoalignment-Finimizers/mutate_seq /Users/biagiele/try/queries/queries_Salmonella_1000_200_${i}.fa /Users/biagiele/try/queries/m_queries_Salmonella_1000_200_${i}_ $m 1

    echo "2. Pseudoalignment with Finimizers"
    # t is not needed now
    #for t in 0.6 0.65 0.7 0.75;
    #do
    #    echo "t=$t"
    #
    #    /usr/bin/time -l /Users/biagiele/try/Pseudoalignment-Finimizers/benchmark search-fmin -i /Users/biagiele/try/Pseudoalignment-Finimizers/results/Salmonella_concat_31 -q /Users/biagiele/try/queries/m_queries_Salmonella_1000_200_${i}_${m}.fa -t $t -o /Users/biagiele/try/Pseudoalignment-Finimizers/stats/m_out_fi_1000_200_${i}_${t}_${m}.txt 1> /Users/biagiele/try/Pseudoalignment-Finimizers/results/search_Salmonella-unitigs-k31.stdout 2> /Users/biagiele/try/Pseudoalignment-Finimizers/results/search_Salmonella-unitigs-k31.stderr
    #done

    #/usr/bin/time -l /Users/biagiele/try/Pseudoalignment-Finimizers/benchmark search-fmin -i /Users/biagiele/try/Pseudoalignment-Finimizers/results/Salmonella_concat_31 -q /Users/biagiele/try/queries/m_queries_Salmonella_1000_200_${i}_${m}.fa -t $t -o /Users/biagiele/try/Pseudoalignment-Finimizers/stats/m_out_fi_1000_200_${i}_${m}.txt 1> /Users/biagiele/try/Pseudoalignment-Finimizers/results/search_Salmonella-unitigs-k31.stdout 2> /Users/biagiele/try/Pseudoalignment-Finimizers/results/search_Salmonella-unitigs-k31.stderr
    /usr/bin/time -l /Users/biagiele/try/Pseudoalignment-Finimizers/benchmark search-fmin -i /Users/biagiele/try/Pseudoalignment-Finimizers/results/Salmonella_concat_31 -q /Users/biagiele/try/queries/queries_Salmonella_1000_200_${i}.fa -t ${t} -o /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_fi_1000_200_${i}.txt 1> /Users/biagiele/try/Pseudoalignment-Finimizers/results/search_Salmonella-unitigs-k31.stdout 2> /Users/biagiele/try/Pseudoalignment-Finimizers/results/search_Salmonella-unitigs-k31.stderr

    #echo "3. Local alignment with BLAST"
    #blastn -db /Users/biagiele/try/Salmonella/db/Salmonella_db -query /Users/biagiele/try/queries/m_queries_Salmonella_1000_200_${i}_${m}.fa -out /Users/biagiele/try/blast/m_out_b_1000_200_${i}_$m.txt
    #rm /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_fi_1000_200_${i}
    #python3 /Users/biagiele/try/parse_blast.py /Users/biagiele/try/blast/m_out_b_1000_200_${i}_${m}.txt /Users/biagiele/try/blast/out/m_out_b_1000_200_${i}_${m}.txt
    # Run Themisto
    #bash /Users/biagiele/try/themisto/run_themisto.sh
    echo "3. Pseudoalignment with Themisto"
    #/Users/biagiele/try/themisto/build/bin/themisto pseudoalign --query-file /Users/biagiele/try/Pseudoalignment-Finimizers/Salmonella/queries/queries_Salmonella_1000_200_${i}.fa --index-prefix /Users/biagiele/try/themisto/my_index --temp-dir /Users/biagiele/try/themisto/temp --out-file /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt --n-threads 1 --threshold 0.7

    /Users/biagiele/try/themisto2/target/release/themisto2 dump-pseudoalignment-data -i /Users/biagiele/try/Salmonella/Salmonella-index.thm2 -q /Users/biagiele/try/queries/queries_Salmonella_1000_200_${i}.fa -m 0 > /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt
    python3 /Users/biagiele/try/parse_json_thm2.py /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt > /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_ti_1000_200_${i}.txt


    # Extract statistics
    echo "4. Extract statistics"
    #for t in 0.6 0.65 0.7 0.75;
    #do
    #    echo "t=$t"

        #python3 /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats.py /Users/biagiele/try/blast/out/m_out_b_1000_200_${i}_$m.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/m_out_fi_1000_200_${i}_${t}_${m}.txt 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats_1000_200_${t}_${m}.txt
        #python3 /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats2.py /Users/biagiele/try/blast/out/m_out_b_1000_200_${i}_$m.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/m_out_fi_1000_200_${i}_${t}_${m}.txt 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats_1000_200_${t}_${m}.txt
    
    #done

    python3 /Users/biagiele/try/Pseudoalignment-Finimizers/compare_rbo.py /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_ti_1000_200_${i}.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_fi_1000_200_${i}.txt 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats_1000_200.txt


    # t = 0.8
    #python3 /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats.py /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_f_1000_200_${i}.stats 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats_1000_200.txt
    #python3 /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats2.py /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_f_1000_200_${i}.stats 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats_1000_200.txt
    
    # intersection
    #python3 /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats.py /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_fi_1000_200_${i}.stats 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/int_stats_1000_200.txt
    #python3 /Users/biagiele/try/Pseudoalignment-Finimizers/stats/stats2.py /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_t_1000_200_${i}.txt /Users/biagiele/try/Pseudoalignment-Finimizers/stats/out_fi_1000_200_${i}.stats 43 >> /Users/biagiele/try/Pseudoalignment-Finimizers/stats/int_stats_1000_200.txt

    echo "DONE $i"
done
echo "Everything done!"
 