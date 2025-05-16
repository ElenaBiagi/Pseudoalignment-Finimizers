.PHONY: finimap
finimap:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++2a -I ./SBWT/sdsl-lite/include/ -O3 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o finimap -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp -lz -D MAX_KMER_LENGTH=250 -D _GLIBCXX_DEBUG -D _GLIBCXX_DEBUG_PEDANTIC

finimap_debug:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o finimap -Wno-deprecated-declarations -march=native -fopenmp -lz -D MAX_KMER_LENGTH=250

tests:
	$(CXX) src/tests.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++2a -I ./SBWT/sdsl-lite/include/ -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o tests -Wno-deprecated-declarations -march=native -fopenmp -lz -D MAX_KMER_LENGTH=250 -export-dynamic
