#include "search_fmin.hh"

int main(int argc, char** argv) {
    string command = argv[1];
    if(command == "search-fmin") {return search_fmin(argc, argv);}
    return 0;
}
// g++ src/test_main.cpp -I ./SBWT/sdsl-lite/include/  -I include -I ./SBWT/include -I ./SBWT/include/sbwt -I SBWT/build/external/sdsl-lite/build/external/libdivsufsort/include/ -o test_main
