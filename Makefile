.PHONY: finimap
SDSL_INCLUDE ?= ./sdsl-lite/include/
SDSL_LIB ?= ./sdsl-lite/build/lib/
finimap:
	$(CXX) src/main.cpp src/utilities.cpp src/seqio.cpp -std=c++2a -O2 -I include -I $(SDSL_INCLUDE) -L $(SDSL_LIB) -march=native -DNDEBUG -fopenmp -lz -lsdsl -D MAX_KMER_LENGTH=250 -o finimap -Wno-deprecated-declarations