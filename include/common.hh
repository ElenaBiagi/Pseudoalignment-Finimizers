#pragma once

#include <string>
#include <cstring>
#include <unordered_map>

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
#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

#include "sbwt/throwing_streams.hh"
#include "PackedStrings.hh"
#include "SeqIO.hh"
//#include "BoundedDeque.hh"




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

// These 3 methods are used in the build phase
pair<int64_t,int64_t> update_sbwt_interval(const int64_t C_char, const pair<int64_t,int64_t>& I, const sdsl::rank_support_v5<>& Bit_rs){
    if(I.first == -1) return I;
    pair<int64_t,int64_t> new_I;
    // both start and end are included
    new_I.first = C_char + Bit_rs(I.first);
    new_I.second = C_char + Bit_rs(I.second+1) -1;
    if(new_I.first > new_I.second){
        return {-1,-1}; // Not found
    } 
    return new_I;
}

pair<int64_t,int64_t> drop_first_char(const int64_t  new_len, const pair<int64_t,int64_t>& I, const sdsl::int_vector<>& LCS, const int64_t n_nodes){
    if(I.first == -1) return I;
    if (new_len<=0){return {0, n_nodes - 1};}
    pair<int64_t,int64_t> new_I = I;
    //Check top and bottom w the LCS
    while (new_I.first > 0 && LCS[new_I.first] >= new_len ){new_I.first --;}
    while(new_I.second < (n_nodes - 1) && LCS[new_I.second + 1] >= new_len ){
        new_I.second ++;
    }
    return {new_I};
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

inline uint64_t suffix2int(const std::string& s, uint64_t offset, char slen) { // if fmin length = 31 we need 62 bits in total, 42 for the tail if plen=10
    uint64_t h = 0;
    for (uint64_t i = 0; i < (uint64_t)slen; i++) {
        uint64_t b = get_char_idx(s[offset + slen - 1 - i]);
        h |= (b << (i << 1)); 
    }
    return h;
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

    for (const auto& [cs, count] : fminFreqCount) {
        //std::cout << cs << " " << count << std::endl; // number of finimizers that appear in x(count) genomes
    }

    for (const auto& [f, count] : fminFreq) {
        //std::cout << f << " " << count << std::endl; // number of genomes with count finimmizers
    }

    
    return;
}

