# Pseudoalignment with Finimizers
Finimizers: [**Finimizers: Variable-length bounded-frequency minimizers for k-mer sets**](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1) by J. N. Alanko, E. Biagi,  S. J. Puglisi. 

### Shortest Unique Finimizers
Let $G$ be the de Bruijn graph of a set of $k$-mers $R$, $t \geq 1$ be an integer, $X$ be a $k$-mer, and $Y$ be a substring of $X$. We say $Y$ is a **shortest $t$-finimizer** of $X$ with respect to the $k$-mer set $R$ if $Y$ has at most $t$ occurrences in $G$ and there does not exist a shorter substring of $X$ with at most $t$ occurrences in $G$. If $t = 1$ then we say $Y$ is a _shortest-unique finimizer_ of $X$.


## Building
First clone the repository and access it with:
```
git clone --recursive https://github.com/ElenaBiagi/Pseudoalignment-Finimizers.git

cd Pseudoalignment-Finimizers
```

If you forgot the flag ```--recursive```, pull the submodules with:
```
git submodule update --init --recursive
```

Then, go to the [finimizer_matrix](https://github.com/jnalanko/finimizer_matrix/tree/2d0127710d8eb6093b43c097de83aaa809da2f6c) submodule and follow the instructions there.
```
cd finimizer_matrix 

cargo build -r

cd ..
```


You are now ready to compile the main project!
```
make finimap
```
## Index construction
Here is a example:

First you should build a colored finimizer matrix. The code takes as input a list of files.

```
Usage: finimizer_matrix build [OPTIONS] --input <INPUT> --output <OUTPUT> --temp-dir <TEMP_DIR> -k <K>
```

```
cd finimizer_matrix 

./target/release/finimizer_matrix build -i ../example_data/coli_file_list.list --reverse -k 31 -t 4 -d ./temp -o ../example_data/coli.cfm 1> /dev/null

```

Then, you can compact the Finimizers index with:
```
build-fmin [OPTION...]

  -i, --index-file arg  ColoredFinimizers file.
  -o, --out-file arg    Output index filename prefix.
  -p, --p_len arg       Finimizers prefix length. (default: 10)
  -x, --x_len arg       Short/long Finimizer threshold. (default: 18)
  -k arg                k-mer length. (default: 31)
  -m, --meta            Metagenome.
  -h, --help            Print usage
```

```
cd ..

./finimap build-fmin -i example_data/coli.cfm  -o example_data/coli_index  -k 31 -p 8 -x 10

```


## Queries

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

```
./finimap search-fmin   -i coli_index.fmin -o example_data/coli_res_t0.txt -q ./example_data/queries.fna -t 0
```

This should output:
```
0 0:50 2:50 1:32 
1 2:50 0:44 1:44 
2 2:26 0:23 1:2
3 2:50 0:21 1:21 
```

## TODO not sorted output at the moment


## Additional info
The code works with the DNA alphabet = {A,C,G,T}.
