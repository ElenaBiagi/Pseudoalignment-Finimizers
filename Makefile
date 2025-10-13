.PHONY: finimap
finimap:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++2a -I ./SBWT/sdsl-lite/include/ -O2 -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -I external/ -I external/bits/external/essentials/include/ -I external/bits/include/ -g -fno-omit-frame-pointer -fno-inline-small-functions -o finimap -Wno-deprecated-declarations -march=native -DNDEBUG -fopenmp -lz -D MAX_KMER_LENGTH=250 -fsanitize=address

finimap_debug:
	$(CXX) src/main.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++20 -I ./SBWT/sdsl-lite/include/ -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o finimap -Wno-deprecated-declarations -march=native -fopenmp -lz -D MAX_KMER_LENGTH=250

tests:
	$(CXX) src/tests.cpp SBWT/build/libsbwt_static.a SBWT/build/external/sdsl-lite/build/lib/libsdsl.a -std=c++2a -I ./SBWT/sdsl-lite/include/ -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -g -o tests -Wno-deprecated-declarations -march=native -fopenmp -lz -D MAX_KMER_LENGTH=250 -export-dynamic
