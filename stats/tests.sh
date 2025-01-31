## Build the index 
# Salmonella   
echo "0.a Build the index for Finimizers"

#/usr/bin/time --verbose ./benchmark build-fmin -u /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/Salmonella_unitigs_31.fna -i /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/results/Salmonella_unitigs_31.sbwt -t 1 -o /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/Salmonella_unitigs_31 -c /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/Salmonella_list.txt 1> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/results/build_Salmonella-unitigs-k31.stdout 2> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/results/build_Salmonella-unitigs-k31.stderr
## --lcs /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/Salmonella_unitigs_31.LCS.sdsl
echo "0.b Build the index for Themisto"
#/home/scratch-hdd/ebiagi/themisto/build/bin/themisto build -k 31 -i /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/Salmonella_list.txt --index-prefix /home/scratch-hdd/ebiagi/themisto/my_index --temp-dir /home/scratch-hdd/ebiagi/themisto/temp --mem-gigas 2 --n-threads 4 --file-colors

rm /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats_100_200.txt
rm /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/int_stats_100_200.txt

for i in $(seq 1 10);
do
    # sample 100 strings of length 200
    echo "1. Sample positive queries"
    #/home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/sample_random /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella_concat.fa 100 200 > /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/queries/queries_Salmonella_100_200_${i}.fa 

    echo "2. Pseudoalignment with Finimizers"

    #/usr/bin/time --verbose /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/benchmark search-fmin -i /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/Salmonella_unitigs_31 -q /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/queries/queries_Salmonella_100_200_${i}.fa -o /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_fi_100_200_${i} 1> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/results/search_Salmonella-unitigs-k31.stdout 2> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/results/search_Salmonella-unitigs-k31.stderr

    #rm /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_fi_100_200_${i}
    # Run Themisto
    #bash /home/scratch-hdd/ebiagi/themisto/run_themisto.sh
    echo "3. Pseudoalignment with Themisto"
    #/home/scratch-hdd/ebiagi/themisto/build/bin/themisto pseudoalign --query-file /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/Salmonella/queries/queries_Salmonella_100_200_${i}.fa --index-prefix /home/scratch-hdd/ebiagi/themisto/my_index --temp-dir /home/scratch-hdd/ebiagi/themisto/temp --out-file /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_t_100_200_${i}.txt --n-threads 1 --threshold 0.7


    # Extract statistics
    echo "4. Extract statistics"
    # t = 0.8
    python3 /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats.py /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_t_100_200_${i}.txt /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_f_100_200_${i}.stats 43 >> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats_100_200.txt
    python3 /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats2.py /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_t_100_200_${i}.txt /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_f_100_200_${i}.stats 43 >> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats_100_200.txt
    
    # intersection
    python3 /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats.py /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_t_100_200_${i}.txt /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_fi_100_200_${i}.stats 43 >> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/int_stats_100_200.txt
    python3 /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/stats2.py /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_t_100_200_${i}.txt /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/out_fi_100_200_${i}.stats 43 >> /home/scratch-hdd/ebiagi/Pseudoalignment-Finimizers/stats/int_stats_100_200.txt

    echo "DONE $i"
done
echo "Everything done!"
 