#include <vector>
#include <iostream>
#include <fstream>
#include "sdsl/bit_vectors.hpp"

using namespace std;

// Colored finimizers without much compression
class ColoredFinimizers {
public:

    vector<char> concat; // Finimizers concatenated in lexicographic order (ASCII characters).
    vector<uint8_t> lengths;
    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors.

    // Loads from the format output by the Rust CLI command `finimizer_matrix` with option --reverse.
    // That format colexicographically sorted reverse finimizers. We reverse them to
    // get lex-sorted finimizers.
    void load(std::istream& in) {
        cerr << "Loading uncompressed tails" << endl;
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

        cerr << "Reversing finimizer strings" << endl;
        int64_t start_in_concat = 0;
        for(int64_t f_idx = 0; f_idx < lengths.size(); f_idx++) {
            int64_t s = start_in_concat;
            int64_t e = start_in_concat + lengths[f_idx];
            std::reverse(concat.begin() + s, concat.begin() + e);
            start_in_concat = e;
        }
    }
};

vector<uint8_t> vbyte_encode(uint64_t x) {
    vector<uint8_t> bytes;
    do {
        uint8_t byte = x & 0x7F;
        x >>= 7;
        if (x != 0) byte |= 0x80;
        bytes.push_back(byte);
    } while (x != 0);
    return bytes;
}

class Bucket {
public:

    // Number of distinct tail lengths in this bucket
    uint8_t n_distinct_lengths; 

    // Compressed tail data
    // Bit layout: [tail length: u8][#tails: vbyte][concat of bitpacked tails]
    //             [tail length: u8][#tails: vbyte][concat of bitpacked tails] 
    //             [tail length: u8][#tails: vbyte][concat of bitpacked tails] 
    //             ...
    sdsl::bit_vector tail_data;

    // Color set ids for each tail
    vector<uint32_t> color_set_ids;

    Bucket(vector<std::string_view> tails, vector<uint32_t>& color_set_ids) {
        // Todo
    }

};

void true_or_crash(bool b, char* error_message){
    if(!b){
        cerr << error_message << endl;
        exit(1);
    }
}

class CompressedColoredFinimizers {

public:

    sdsl::bit_vector color_sets_concat; // Length #finimizers * #colors. TODO: deduplicate.
    vector<Bucket> buckets;
    uint64_t n_colors;
    uint64_t n_finimizers;

    CompressedColoredFinimizers(ColoredFinimizers& cf, int64_t prefix_len) {
        uint64_t n_finimizers = cf.lengths.size();
        true_or_crash(n_finimizers > 0, "ERROR: 0 finimizers");

        vector<std::string_view> cur_bucket_nonnegative_tails;
        int64_t first_nonegative_tail_idx = -1;
        int64_t f_start = 0;
        for(int64_t i = 0; i < n_finimizers; i++){
            if(cf.lengths[i] >= prefix_len) {
                first_nonegative_tail_idx = i;
                break;
            }
            f_start += cf.lengths[i];
        }
        true_or_crash(first_nonegative_tail_idx >= 0, "ERROR: all tails shorter than prefix length");

        std::string_view cur_prefix(cf.concat.data() + f_start, prefix_len);
        vector<std::string_view> cur_tails;
        vector<uint32_t> cur_color_set_ids;
        for(int64_t i = 0; i < n_finimizers; i++) {
            if(cf.lengths[i] < prefix_len){
                // Skip negative tails. TODO: do something about them.
            } else {
                std::string_view prefix(cf.concat.data() + f_start, prefix_len);
                if(prefix != cur_prefix) {
                    // Bucket changes -> encode currently collected tails
                    buckets.push_back(Bucket(cur_tails, cur_color_set_ids));
                    cur_tails.clear();
                    cur_color_set_ids.clear();
                }
                cur_tails.push_back(std::string_view(cf.concat.data() + f_start + prefix_len, cf.lengths[i] - prefix_len));
                cur_color_set_ids.push_back(i);
                cur_prefix = prefix;
            }
            f_start += cf.lengths[i];
        }

        if(cur_tails.size() > 0){ // Last bucket
            buckets.push_back(Bucket(cur_tails, cur_color_set_ids));
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
