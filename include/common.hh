#pragma once

#include <string>
#include <cstring>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
#include <tuple>
#include <set>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>


#include "sbwt/cxxopts.hpp"
#include "sbwt/globals.hh"
#include "sbwt/SBWT.hh"
#include "sbwt/SubsetWT.hh"
#include "sbwt/stdlib_printing.hh"
#include "sbwt/SeqIO.hh"
#include "sbwt/SubsetMatrixRank.hh"
#include "sbwt/buffered_streams.hh"
#include "sbwt/variants.hh"
#include "sbwt/commands.hh" 

#include "SeqIO.hh"

using namespace std;
using namespace sbwt;

/* // Table mapping ascii values of characters to their reverse complements,
// lower-case to lower case, upper-case to upper-case. Non-ACGT characters
// are mapped to themselves.
static constexpr unsigned char rc_table[256] =
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73,
74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91,
92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122,
123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242,
243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};


// ACGT -> TGCA
constexpr char get_rc(char c){
    return rc_table[(unsigned char)c];
}
 */
set<string> get_substrings(const string& str) {
    set<string> substrings;
    int n = str.size();
    for (int i = 0; i < n; ++i) {
        for (int len = 1; len <= n - i; ++len) {
            substrings.insert(str.substr(i, len));
        }
    }
    return substrings;
}

inline void print_bit_vector(const std::unique_ptr<sdsl::bit_vector>& T) {
    if (!T) {
        std::cout << "T is null." << std::endl;
        return;
    }

    std::cout << "bit_vector of size " << T->size() << ": [ ";
    for (size_t i = 0; i < T->size(); ++i) {
        std::cout << (*T)[i] << " ";
        if (i == 128) { break; } 
    }
    std::cout << "]" << std::endl;
}


// TODO remove
inline void print_results(const std::unordered_map<int, uint64_t>& results) {
    std::cout << "results (position → count):" << endl;
    for (const auto& [pos, count] : results) {
        std::cout << "  " << pos << " → " << count << endl;
    }
}

inline char get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// TODO remove
uint64_t prefix2int_old(const string& s, uint64_t offset, char plen){ 
    uint64_t h = 0;
    for(uint64_t i=0; i<(uint64_t)plen; i++){
       uint64_t b = get_char_idx(s[i+offset]);
       h |= (b << (i<<1));
    }
    //cerr << h << '\n';
    return h;
}

inline uint64_t prefix2int(const std::string& s, uint64_t offset, char plen){ // if fmin length = 31 we need 62 bits in total, 20 for the prefix if plen=10
    uint64_t h = 0;
    for (uint64_t i = 0; i < (uint64_t)plen; i++) {
        uint64_t b = get_char_idx(s[offset + i]);
        h <<= 2;
        h |= b;
    }
    return h;
}

inline uint64_t prefix2int(const std::string_view s, uint64_t offset, char plen) {
    uint64_t h = 0;
    for (uint64_t i = 0; i < static_cast<uint64_t>(plen); ++i) {
        uint64_t b = get_char_idx(s[offset + i]);
        h <<= 2;
        h |= b;
    }
    return h;
}

inline uint64_t stream_kmer(uint64_t prev_hash, char new_char, char plen) {
    uint64_t b = get_char_idx(new_char);
    prev_hash <<= 2;              // shift by 2
    prev_hash |= b;               // new char
    prev_hash &= ((1ULL << (2 * plen)) - 1); // keep only plen bases (mask older bits)
    return prev_hash;
}


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
    
inline sdsl::int_vector<1> WriteTailsVector(vector<vector<string>>& tails,vector <int> tlens ){
    uint64_t total_bits = 0;
    for (size_t i = 0; i < tlens.size(); ++i) {
        int tlen = tlens[i];
        uint32_t ntails = tails[i].size();

        total_bits += 5; // tlen
        total_bits += vbyte_encode(ntails).size() * 8; // vbyte encoded tail count

        total_bits += ntails * tlen * 2; // each tail uses tlen*2 bits
    }

    sdsl::int_vector<1> T;
    T.resize(total_bits+128); 
    uint64_t* data = T.data();


   const char ntlen = tlens.size();
    int64_t offset = 0;
    uint64_t word_index = 0;
    uint8_t w_offset = 0;
           
    for (char i=0; i<ntlen; i++){
        int tlen = tlens[i]; 
        word_index = offset/64;
        w_offset = offset %64;
        sdsl::bits::write_int(&data[word_index], tlen, w_offset, 5);
        offset += 5;

        const uint32_t tnumber = tails[i].size();
        auto vb = vbyte_encode(tnumber);
        for (uint8_t b : vb) {
            word_index = offset/64;
            w_offset = offset %64;
            sdsl::bits::write_int(&data[word_index], b, w_offset, 8);
            offset += 8;
        }
        for (const auto &t : tails[i]){ //pair<uint32_t, set<int>
            // tails
            word_index = offset/64;
            w_offset = offset %64;
            uint32_t t_int = prefix2int(t,0,tlen);
            sdsl::bits::write_int(&data[word_index], t_int, w_offset, tlen*2);
            offset += (2 * tlen);
        }
    }
    return T;
}


// TODO REMOVE
//used in build-verify
vector<string> remove_ns(const string& unitig, const int64_t k){
    vector<string> new_unitigs;
    const int64_t str_len = unitig.size();
    int64_t start = 0;
    char c;
    char char_idx;
    for (int64_t i = 0; i < str_len;i++){
        c = static_cast<char>(unitig[i] &~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]] {
            if ((i - start + 1) >= k ){
                string new_seq = unitig.substr(start,(i - start + 1));
                new_unitigs.push_back(new_seq);
            }
            start = i + 1;
        }
    }
    if ((str_len - start) >= k ){
        string new_seq = unitig.substr(start,(str_len - start));
        new_unitigs.push_back(new_seq);
    }
    return new_unitigs;
}

//not used
// TODO REMOVE
const std::string remove_N_from_string(const std::string &input) {
    std::string result = input;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper); // uppercase
    result.erase(std::remove(result.begin(), result.end(), 'N'), result.end());
    return result; // Return the resulting string as const
}

// TODO REMOVE
vector< std::string> split_by_N(const std::string &input, const int64_t k) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t end = 0;

    while ((end = input.find_first_of("BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrstuvwxyz", start)) != std::string::npos) {
        if (end - start +1 >= k) { // Exclude strings shorter than k
            string seq = input.substr(start, end - start);
            std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
            result.emplace_back(seq);
        }
        start = end + 1;
    }

    // Add the last part if it's at least k characters long
    if (input.size()-start >=k) {
        string seq = input.substr(start);
        std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        result.emplace_back(seq);
    }

    return result;
}

string print_finimizer_stats(const vector<tuple<int64_t, int64_t, int64_t>>& finimizers, int64_t n_kmers, int64_t t, const sdsl::bit_vector& color_sets_concat){
    // len, int, colors
    int64_t new_number_of_fmin = finimizers.size();
    set<tuple<int64_t, int64_t, int64_t>> set_fmin = set(finimizers.begin(), finimizers.end());
    int64_t d_number_of_fmin = set_fmin.size();
    
    vector<uint64_t> lengths;
    lengths.resize(31);

    uint64_t n_colors = 1990;
    uint64_t sum_colors = 0;
    uint64_t s_colors = 0;

    vector<uint64_t> v_colors;
    v_colors.resize(1990);

    //int64_t sum_freq = 0;

    int64_t sum_len = 0;
    for (auto x : finimizers){
        s_colors = 0;
        //sum_freq += get<1>(x);
        auto start = get<2>(x);
        for (uint64_t i=0; i< n_colors; i++){
            s_colors+=color_sets_concat[(n_colors*start)+i];
        }
        sum_colors += s_colors;
        v_colors[s_colors] += 1;
        lengths[get<0>(x)] += 1;
        sum_len += get<0>(x);
    }

    string result = to_string(new_number_of_fmin) + "," + to_string(d_number_of_fmin) + "," + to_string(static_cast<float>(sum_colors) / static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len) / static_cast<float>(new_number_of_fmin));

    write_log(to_string(t) + "," + result, LogLevel::MAJOR);
    write_log("#total finimizers: " + to_string(new_number_of_fmin) , LogLevel::MAJOR);

    write_log("#Distinct finimizers: " + to_string(d_number_of_fmin) , LogLevel::MAJOR);

    //write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
    //write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    write_log("Avg colors: " + to_string(static_cast<float>(sum_colors)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);

    write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    return result;
}

    std::ostream& operator<<(std::ostream& os, const std::set<int>& set) {
    os << "{";
    for (auto it = set.begin(); it != set.end(); ++it) {
        os << *it;
        if (std::next(it) != set.end()) {
            os << ", ";
        }
    }
    os << "}";
    return os;
}
