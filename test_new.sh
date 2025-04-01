folder=$1
colors=$2

if [ ! -d "$folder" ]; then
  echo "$folder does not exist."
  exit
fi
echo $folder

echo "0. Build all"
./build_all.sh $folder


echo "1. Sample pos queries"
./extract_queries.sh $folder

echo "3. Query all"
./query_all.sh $folder

echo "4. Extract statistics"
./stats_all.sh $folder $colors
