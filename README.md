# Pseudoalignment with Finimizers
Finimizers: [**Finimizers: Variable-length bounded-frequency minimizers for k-mer sets**](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1) by J. N. Alanko, E. Biagi,  S. J. Puglisi. 

### Shortest Unique Finimizers
Let $G$ be the de Bruijn graph of a set of $k$-mers $R$, $t \geq 1$ be an integer, $X$ be a $k$-mer, and $Y$ be a substring of $X$. We say $Y$ is a **shortest $t$-finimizer** of $X$ with respect to the $k$-mer set $R$ if $Y$ has at most $t$ occurrences in $G$ and there does not exist a shorter substring of $X$ with at most $t$ occurrences in $G$. If $t = 1$ then we say $Y$ is a _shortest-unique finimizer_ of $X$.


## Building
First clone the repository with:
```
git clone --recursive https://github.com/ElenaBiagi/Pseudoalignment-Finimizers.git 
```

If you forgot the flag ```--recursive```, pull the submodules with:
```
git submodule update --init --recursive
```

Then, go the [SBWT](https://github.com/algbio/SBWT/tree/eb7f54165d38bb7c9aa2ab418b7d17a1113a9977) submodule and build it using the instructions in the submodule. And compile the experiments with:
```
cd SBWT/build

cmake .. -DCMAKE_C_COMPILER=$(which gcc-10) -DCMAKE_CXX_COMPILER=$(which g++-10) -D MAX_KMER_LENGTH=250
make -j4

cd ../..
```
Then, go to the [finimizer_matrix](https://github.com/jnalanko/finimizer_matrix/tree/2d0127710d8eb6093b43c097de83aaa809da2f6c) submodule and follow the instructions there.


You are now ready to compile the main project!
```
make benchmark --always-make CXX=g++-10
```
## Index construction


First ypu should build a finimizer matrix. The code takes as input a list of files.
Here is a example:
# Fix example

```
cd Pseudoalignment-Finimizers/finimizer_matrix
./target/release/finimizer_matrix build -i <files_list.txt> --reverse -k 31 -t 12 -o Salmonella.cfm -d ./temp -m 250

```

Then, you can compact the Finimizers index with:

```
./finimap build-fmin -o <colored-finimizer-index>  -i Salmonella.cfm -k 31 -p 10
```
```
Usage:
  build-fmin [OPTION...]

  -i, --index-file arg  ColloredFinimizers file.
  -o, --out-file arg    Output index filename prefix.
  -p, --p_len arg       Finimizers prefix length. (default: 10)
  -k arg                k-mer length. (default: 31)
  -h, --help            Print usage
```

## Queries

```
./finimap search-fmin -o <out-file>  -i <colored-finimizer-index> -q <query-file.fa> -t <threshold>
```
```
Usage:
  search-fmin [OPTION...]

  -o, --out-file arg    Output filename, or stdout if not given.
  -i, --index-file arg  Index filename prefix.
  -q, --query-file arg  The query in FASTA or FASTQ format, possibly 
                        gzipped. Multi-line FASTQ is not supported. If the 
                        file extension is .txt, this is interpreted as a 
                        list of query files, one per line. In this case, 
                        --out-file is also interpreted as a list of output 
                        files in the same manner, one line for each input 
                        file.
  -t arg                Threshold (default: 0)
  -h, --help            Print usage
```
The result of a query will be the number or percentage (t > 0) of finimizers observed per color, expressed in pairs of (color:#matches).

## RBO


## Additional info
The code works with the DNA alphabet = {A,C,G,T}.
