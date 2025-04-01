# Remember to give in input the number of colors
folder=$1
colors=$2

if [ ! -d "$folder/stats" ]; then
  mkdir ${folder}/stats
fi

m=0.01
for i in $(seq 1 10);
do  
    # Themisto vs Finimap
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_ti_1000_200_${i}.txt $HOME/${folder}/results/out_fi_1000_200_${i}.txt 43 >> $HOME/${folder}/stats/stats_ft_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_ti_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_fi_1000_200_${i}_${m}.txt 43 >> $HOME/${folder}/stats/m_stats_ft_1000_200_${m}.txt
    
    # Themisto vs Themisto
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_ti_1000_200_${i}.txt $HOME/${folder}/results/out_ti_1000_200_${i}.txt $colors >> $HOME/${folder}/stats/stats_tt_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_ti_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_ti_1000_200_${i}_${m}.txt $colors >> $HOME/${folder}/stats/m_stats_tt_1000_200_${m}.txt

    # COBS vs Finimap
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_ci_1000_200_${i}.txt $HOME/${folder}/results/out_fi_1000_200_${i}.txt 43 >> $HOME/${folder}/stats/stats_fc_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_ci_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_fi_1000_200_${i}_${m}.txt 43 >> $HOME/${folder}/stats/m_stats_fc_1000_200_${m}.txt
    
    # BLAST vs Finimap
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_bi_1000_200_${i}.txt $HOME/${folder}/results/out_fi_1000_200_${i}.txt 43 >> $HOME/${folder}/stats/stats_fb_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_bi_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_fi_1000_200_${i}_${m}.txt 43 >> $HOME/${folder}/stats/m_stats_fb_1000_200_${m}.txt
    
    # Themisto vs COBS
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_ti_1000_200_${i}.txt $HOME/${folder}/results/out_ci_1000_200_${i}.txt $colors >> $HOME/${folder}/stats/stats_tc_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_ti_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_ci_1000_200_${i}_${m}.txt $colors >> $HOME/${folder}/stats/m_stats_tc_1000_200_${m}.txt

    # Themisto vs BLAST
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_ti_1000_200_${i}.txt $HOME/${folder}/results/out_bi_1000_200_${i}.txt $colors >> $HOME/${folder}/stats/stats_tb_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_ti_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_bi_1000_200_${i}_${m}.txt $colors >> $HOME/${folder}/stats/m_stats_tb_1000_200_${m}.txt

    # COBS vs BLAST
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/out_ci_1000_200_${i}.txt $HOME/${folder}/results/out_bi_1000_200_${i}.txt $colors >> $HOME/${folder}/stats/stats_cb_1000_200.txt
    python3 $HOME/Pseudoalignment-Finimizers/compare_rbo.py $HOME/${folder}/results/m_out_ci_1000_200_${i}_${m}.txt $HOME/${folder}/results/m_out_bi_1000_200_${i}_${m}.txt $colors >> $HOME/${folder}/stats/m_stats_cb_1000_200_${m}.txt
done
