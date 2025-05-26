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

inline void print_bit_vector(const sdsl::bit_vector& T, int64_t pos) {
    std::cerr << "T size = " << T.size() << std::endl;
    std::cerr << "pos = " << pos << std::endl;
    if (pos + 4 < (int64_t)T.size()) {
        std::cerr << "T[pos] = " << T[pos] << " " << T[pos+1] << " " << T[pos+2] << " " << T[pos+3] << " " << T[pos+4] << std::endl;
    } else {
        std::cerr << "pos too close to end of vector for T[pos+4]" << std::endl;
    }

    std::cout << "bit_vector content from pos: ";
    for (size_t i = pos; i < T.size() && i < (size_t)(pos + 64); ++i) {
        std::cout << T[i] << " ";
    }
    std::cout << std::endl;
}

void print_B(const unordered_map<uint32_t, pair<int64_t,int64_t> >& B) {
    std::cout << "B (prefix → offset): " << endl;
    for (const auto& [prefix_hash, offset] : B) {
        std::cout << prefix_hash << " → { " << offset.first << ", " << offset.second << " }" << endl;
    }
}

void print_helperB(const std::unordered_map<uint32_t, std::map<char, std::set<std::pair<uint32_t, std::set<int>>>>>& helperB) {
    std::cerr << "helperB:" << endl;
    for (const auto& [prefix, tail_map] : helperB) {
        std::cerr << "Prefix: " << prefix << endl;
        for (const auto& [character, finimizer_set] : tail_map) {
            std::cerr << "  └─ Tlen: '" << (int)character << "' → " << finimizer_set.size() << " finimizer(s)\n";
            for (const auto& [tail, color_set] : finimizer_set) {
                std::cerr << "      └─ Tail: " << tail << " → Colors: { ";
                for (int color : color_set) {
                    std::cerr << color << " ";
                }
                std::cerr << "}"<<endl;
            }
        }
    }
}

void print_sB(const std::unordered_map<uint32_t, int64_t>& sB) {
    std::cout << "sB (fmin → offset):\n";
    for (const auto& [key, value] : sB) {
        std::cout << key << " → " << value << endl;
    }
}

void printHashTable(const std::unordered_map<std::string, std::set<int>>& hashTable) {
    std::cerr << "HASH TABLE" << std::endl;
    for (const auto& pair : hashTable) {
        std::cerr << "Key: " << pair.first << ", Values: " << pair.second << std::endl;
    }
    std::cerr << "HASH TABLE done" << std::endl;
}

void print_results(const std::unordered_map<int, uint64_t>& results) {
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

inline uint64_t suffix2int(const std::string& s, uint64_t offset, char slen) { // if fmin length = 31 we need 62 bits in total, 42 for the tail if plen=10
    uint64_t h = 0;
    for (uint64_t i = 0; i < (uint64_t)slen; i++) {
        uint64_t b = get_char_idx(s[offset + slen - 1 - i]);
        h |= (b << (i << 1)); 
    }
    return h;
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


//old
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
/* void remove_N_from_string(std::string &s) {
    s.erase(std::remove(s.begin(), s.end(), 'N'), s.end());
} */
const std::string remove_N_from_string(const std::string &input) {
    std::string result = input;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper); // uppercase
    result.erase(std::remove(result.begin(), result.end(), 'N'), result.end());
    return result; // Return the resulting string as const
}

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

//TODO we might want to print the stats of finimizers found in all the genomes
void get_stats(std::unordered_map<std::string, std::set<int>>& hashTable){  
    //std::unordered_map<std::string, int> genomes;
    std::map<int, int> fminFreqCount;      
    std::map<int, int> fminFreq;

    std::map<int, int> freq;  // number of finimizers for each genome
    std::cout << hashTable.size()<< std::endl;
    for (const auto& [fmin, colors] : hashTable) {
        //genomes[fmin] = colors.size();
        size_t cs = colors.size();
        fminFreqCount[cs]++;
        //std::cout << fmin << ": " << cs << std::endl;
        for (int c : colors) { // add +1 to every color/genome observed 
            freq[c]++;
        }
    }
    
    for (const auto& [c, f] : freq) {
        fminFreq[f]++; // number of genomes with f finimizers
        std::cout << c << ": " << f << std::endl;
    }

/*     for (const auto& [cs, count] : fminFreqCount) {
        //std::cout << cs << " " << count << std::endl; // number of finimizers that appear in x(count) genomes
    }

    for (const auto& [f, count] : fminFreq) {
        //std::cout << f << " " << count << std::endl; // number of genomes with count finimmizers
    } */
    return;
}