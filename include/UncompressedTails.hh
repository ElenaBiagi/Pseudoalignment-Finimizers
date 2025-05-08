#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "sdsl/bit_vectors.hpp"

using namespace std;

class UncompressedTails {
public:

    vector<char> concat; // Ascii characters
    vector<uint8_t> lengths;
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors

    void load(std::istream& in) {
        uint64_t n_finimizers;
        in.read(reinterpret_cast<char*>(&n_finimizers), sizeof(n_finimizers));

        uint64_t finimizer_total_length;
        in.read(reinterpret_cast<char*>(&finimizer_total_length), sizeof(finimizer_total_length));

        uint64_t n_colors;
        in.read(reinterpret_cast<char*>(&n_colors), sizeof(n_colors));

        lengths.resize(n_finimizers);
        in.read(reinterpret_cast<char*>(lengths.data()), n_finimizers * sizeof(uint8_t));

        concat.resize(finimizer_total_length);
        in.read(reinterpret_cast<char*>(concat.data()), finimizer_total_length * sizeof(char));

        int64_t n_bits = n_finimizers * n_colors;
        // The bits are in u64 Lsb format
        vector<uint64_t> words((n_bits + 63) / 64); // Ceil div by 64
        in.read(reinterpret_cast<char*>(words.data()), words.size() * sizeof(uint64_t));
        color_sets_concat.resize(n_bits);
        for(int64_t i = 0; i < words.size(); i++) {
            color_sets_concat.set_int(i*64, words[i]);
        }
    }
};

/*
int main(){
    UncompressedTails ut;
    ifstream in("out.bin");
    ut.load(in);
    for(auto c : ut.concat) cout << c; cout << endl;
    for(auto len : ut.lengths) cout << (int)len << " "; cout << endl;
    cout << ut.color_sets_concat << endl;
}
*/