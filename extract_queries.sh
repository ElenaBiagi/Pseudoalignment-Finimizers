if [ ! -d "$folder/queries" ]; then
    mkdir ${folder}/queries
fi  

for i in $(seq 1 10);
do
    echo "1. Sample positive queries ($i)"
    $HOME/Pseudoalignment-Finimizers/sample_random $HOME/${folder}/${folder}_concat.fa 1000 200 i > $HOME/${folder}/queries/queries_${folder}_1000_200_${i}.fa 

    echo "2. Mutate queries"
    $HOME/Pseudoalignment-Finimizers/mutate_seq $HOME/queries/queries_${folder}_1000_200_${i}.fa $HOME/${folder}/queries/m_queries_${folder}_1000_200_${i}_ $m i
done