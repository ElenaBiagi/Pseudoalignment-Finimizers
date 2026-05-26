# Finimap

**Finimap** is a tool for building and querying compressed colored finimizer indexes for fast sequence pseudoalignment.

Original paper: **Finimap: fast and accurate single-species bacterial pseudoalignment with finimizers** by J. N. Alanko, E. Biagi, S. J. Puglisi.

Finimizers paper: [**Finimizers: Variable-length bounded-frequency minimizers for k-mer sets**](https://www.biorxiv.org/content/10.1101/2024.02.19.580943v1) by J. N. Alanko, E. Biagi, S. J. Puglisi.

## Overview

A finimizer index enables fast querying of sequence data (e.g., DNA reads) against reference genomes or genome collections. The tool builds a compact index by:

1. **Extracting finimizers** from reference sequences
2. **Compacting the index**
3. **Querying efficiently** by extracting finimizers from query sequences and finding matches in the index

This approach provides superior speed compared to traditional k-mer indexes while maintaining query accuracy.

### Shortest Unique Finimizers

Let $G$ be the de Bruijn graph of a set of $k$-mers $R$, $t \geq 1$ be an integer, $X$ be a $k$-mer, and $Y$ be a substring of $X$. We say $Y$ is a **shortest $t$-finimizer** of $X$ with respect to the $k$-mer set $R$ if $Y$ has at most $t$ occurrences in $G$ and there does not exist a shorter substring of $X$ with at most $t$ occurrences in $G$. If $t = 1$ then we say $Y$ is a _shortest-unique finimizer_ of $X$.


## Building

### Prerequisites

- Rust compiler (for building the `finimizer_matrix` submodule)
- A C/C++ compiler (for the main project)
- Make

### Compilation Steps

First clone the repository and access it with:
```bash
git clone --recursive https://github.com/ElenaBiagi/Pseudoalignment-Finimizers.git
cd Pseudoalignment-Finimizers
```

If you forgot the flag `--recursive`, pull the submodules with:
```bash
git submodule update --init --recursive
```

Build the [finimizer_matrix](https://github.com/jnalanko/finimizer_matrix/tree/2d0127710d8eb6093b43c097de83aaa809da2f6c) submodule:
```bash
cd finimizer_matrix
cargo build -r
cd ..
```

Build the [sdsl-lite](https://github.com/iosfwd/sdsl-lite.git) dependency:
```bash
cd sdsl-lite && mkdir -p build && cd build
cmake .. -DCMAKE_CXX_FLAGS="-std=c++2a -O2 -march=native"
make -j4
cd ../..
```

Compile the main project:
```bash
make finimap
```

This will produce two main executables:
- `finimizer_matrix` (in the `finimizer_matrix/target/release/` directory)
- `finimap` (in the current directory)
## Index Construction

Building a finimizer index is a two-stage process:

1. **Stage 1**: Build a colored finimizer matrix (CFM) from reference sequences
2. **Stage 2**: Compact the CFM into a queryable finimizer index (FMIN)

### Stage 1: Building the Colored Finimizer Matrix

The first stage extracts finimizers from your reference sequences using `finimizer_matrix build`.

#### Input Format

Create a list file (`.list`) containing paths to your reference sequence files (one per line):

```bash
# example_data/coli_file_list.list
../example_data/GCA_000005845.2_ASM584v2.fna
../example_data/GCA_000006665.1_ASM666v1.fna
../example_data/GCA_000007445.1_ASM744v1.fna
```

#### Command

```bash
cd finimizer_matrix

./target/release/finimizer_matrix build [OPTIONS] --input <INPUT_LIST> --output <OUTPUT_CFM> --temp-dir <TEMP_DIR> -k <K_VALUE>
```

#### Parameters

- `-i, --input <FILE>`: Path to list file containing reference sequence paths (one per line, FASTA format)
- `-o, --output <FILE>`: Output colored finimizer matrix file
- `-k <K>`: K-mer length (default: 31). Recommended values: 19-31
- `-t <THREADS>`: Number of threads (default: 1)
- `-d, --temp-dir <DIR>`: Temporary directory for intermediate files
- `-m <MEMORY>`: Maximum memory usage in MB
- `--reverse`: Include reverse complement of sequences 
- `--sparse`: Use sparse representation for finimizers (saves memory/space)

#### Example

```bash
cd finimizer_matrix

./target/release/finimizer_matrix build \
  -i ../example_data/coli_file_list.list \
  -o ../example_data/coli.cfm \
  -k 31 \
  -t 4 \
  -d ./temp \
  --reverse \
  2> /tmp/build.log
```

This builds a CFM file named `coli.cfm` from the references listed in `coli_file_list.list`.

### Stage 2: Compacting the Index with `build-fmin`

The second stage compacts the CFM into an efficient queryable index format.

#### Command

```bash
cd ..

./finimap build-fmin -i <CFM_FILE> -o <INDEX_PREFIX> [OPTIONS]
```

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-i, --index-file` | FILE | Required | Colored finimizer matrix (CFM) file from Stage 1 |
| `-o, --out-file` | FILE | Required | Output index filename prefix (without extension) |
| `-p, --p_len` | INT | 10 | Finimizer prefix length (bits used for bucketing). Higher values → more buckets → faster but larger index |
| `-x, --x_len` | INT | 18 | Short/long finimizer threshold (bits). Finimizers shorter than this use one representation, longer use another |
| `-k` | INT | 31 | K-mer length (must match the value used in Stage 1) |
| `-m, --meta` | BOOL | false | Set to `true` for metagenomic indexes (multiple reference genomes) |
| `-h, --help` | - | - | Display help message |

#### Parameter Tuning

- **`-p` (prefix length)**: Controls index bucket size in long finimizers. Typical values: 8-12, (p > x)
  - Smaller values: smaller index, slower queries
  - Larger values: faster queries, larger index size
  
- **`-x` (short/long threshold)**: Affects short/long finimizer representation
  - Typical values: 10-20
  - Adjust based on your finimizer distribution

#### Examples

**Single genome index (e.g., *E. coli*):**
```bash
./finimap build-fmin \
  -i example_data/coli.cfm \
  -o example_data/coli_index \
  -k 31 \
  -p 10 \
  -x 18
```

This produces an index `example_data/coli_index.fmin`.

**Compact index with smaller memory footprint:**
```bash
./finimap build-fmin \
  -i example_data/coli.cfm \
  -o example_data/compact_index \
  -k 31 \
  -p 8 \
  -x 15
```


## Querying

The `search-fmin` command finds finimizer matches in your query sequences against the built index.

### Command

```bash
./finimap search-fmin -i <INDEX_PREFIX> -q <QUERY_FILE> -o <OUTPUT_FILE> [OPTIONS]
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `-i, --index-file` | FILE | Required | Index filename prefix (from `build-fmin`) |
| `-q, --query-file` | FILE | Required | Query sequences in FASTA/FASTQ format (may be gzipped). Or `.txt` file with list of query files (one per line) |
| `-o, --out-file` | FILE | Required | Output filename (or stdout if omitted). Or `.txt` file with output file list (one per line, paired with query list) |
| `-t` | INT | 0 | Threshold for reporting matches in threshold union queries. If 0: absolute counts; if > 0: percentage of matching k-mers; if 1 intersection union |
| `-b` | INT | 1000 | Number of batched queries. |
| `-h, --help` | - | - | Display help message |

### Output Format

The output contains one line per query sequence:
```
QUERY_ID color_1:count_1 color_2:count_2 color_3:count_3 ...
```

Where:
- **QUERY_ID**: Sequence identifier from the query file (int)
- **color_N**: Index of a reference genome in the index
- **count_N**: Number of matching finimizers (absolute or percentage depending on `-t` threshold)

Colors are assigned in the order references were listed in the original input file.

```color:count``` pairs are sorted in decreasing order. 



### Examples

**Threshold union queries (no minimium match requirement):**
```bash
./finimap search-fmin \
  -i example_data/coli_index \
  -q example_data/queries.fna \
  -b 4 \
  -o results.txt \
  -t 1
```

Output example:
```
0 0:50 2:50 1:32 
1 2:50 0:44 1:44 
2 2:26 0:23 1:2 
3 2:50 0:21 1:21
```

**Threshold union queries (80% match required):**
```bash
./finimap search-fmin \
  -i example_data/coli_index \
  -q example_data/queries.fna \
  -b 4 \
  -o results.txt \
  -t 0.8
```


## Additional Information

- **Alphabet**: The code works with DNA sequences: {A, C, G, T} (case-insensitive)
- **Reverse complement**: Use `--reverse` flag in `finimizer_matrix build` to index both forward and reverse complement strands
- **Metagenomic data**: Use `--meta true` in `build-fmin` for collections of multiple genomes
- **Large databases**: For very large reference collections (>100GB), consider:
  - Using higher `-p` values for faster queries
  - Using the sparse representation