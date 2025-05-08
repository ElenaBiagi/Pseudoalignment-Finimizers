#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

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
#include "BoundedDeque.hh"

#include "common.hh"
#include "bitsearch.hh"


// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? YES
// set ?
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const unordered_map<uint32_t, pair<int64_t,int64_t> >& B, const std::unordered_map<uint32_t, int64_t>& sB, const std::vector<int_vector<1>>& T, const uint8_t plen, const int k){ 
    
    const int64_t str_len = input.size();

    vector<uint64_t> Fmin;// pointer to C

    int64_t last_pos = 0;
    int64_t start = 0;
    int64_t kmer_start = 0; // start of the first k-mer

    bool found = false;
    uint8_t len_fmin = 0;
    uint64_t int_fmin = 0;
    //string str_fmin;
    uint64_t s_int = 0;
    int64_t C_fmin = 0; // Can the result of color index be negative??? If not found??


    BoundedDeque<tuple<uint8_t, uint64_t, int64_t, int64_t>> all_fmin(input.size()-k+1);
    uint64_t m = std::numeric_limits<uint64_t>::max();

    tuple<uint8_t, uint64_t, int64_t, int64_t> curr_substr; // length, fmin, C_offset, start
    tuple<uint8_t, uint64_t, int64_t, int64_t> w_fmin = {k+1,m,0, input.size()}; // start will always be < str_len
    
    // idea: look for prefixes of length p in the hashtable B
    // 1. prefix not found: the finimizer might be smaller
    //     Check sB with a shorter prefix
    //      * found = store finimizer
    //      - not found = continue
    // 2. prefix found
    //    Go to where the pointer takes you in T
    //    Check first the shortest lengths

        
    // TODO select the correct fmin for every k-mer
    for (start = 0; start < str_len - k +1; start++) { // Extract p characters at a time
        // TODO convert 32 values at time = 64 bits
        uint64_t int_p = prefix2int(input, start, plen);
        // Look for the prefix in B
        //if (B.find(intp) != B.end()) { // B contains all the possible prefixes of length p
        
        //int64_t pointer = B.find((uint32_t)int_p)->second;
        auto pp = B.at(int_p);
        int64_t pointer = pp.first;
        int64_t tails_so_far = pp.second;


        if (pointer != -1){
            // 2. Prefix found!
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch_new depending on tlen
            char s_len = (str_len >= start+k) ? k-plen : str_len-start;
            s_int = prefix2int(input, start+plen, s_len);
            auto result = bitMagicSearch_new(T[pointer], s_int, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            
            if (result.first != -1){ 
                // b. Tail Found!
                found = true;
                len_fmin = plen + result.second; // TODO add the correct fmin len
                C_fmin = result.first + tails_so_far; // TODO NO NEED TO STORE THE COLORS NOW AS LONG AS WE KEEP THE OFFSET 
                // Store the string as a number. OK as only strings of the same length will be compared 
                
                // TODO extract the number instead of converting again
                int_fmin = (s_int >> ((s_len - len_fmin) * 2)) & ((1ULL << (len_fmin * 2)) - 1); 
                //int_fmin =  prefix2int(input, start+plen, len_fmin);
                
                // TODO compare LEXICOGRAPHICALLY
                
            }
        } else{
            // 1. Prefix NOT found
            // Shorten the prefix until you find a match
            // TODO should we keep the length of the SHORTEST finimizer? To know when to stop
            uint8_t sp_len = plen-1;
            //uint64_t int_sp = prefix2int(input, start, slen); // It might be smarter to start from the longest prefix, done
            uint64_t int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1); // subtract 2 bits from the original prefix
            auto it = sB.find(int_sp);

            while(it == sB.end() and sp_len > 0){
                sp_len--;
                int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1); // Subtract 2bits (1 letter) at a time
                //int_sp = prefix2int(input, start, sp_len);
            }
            if (sp_len != 0){
                found = true;
                len_fmin = (uint8_t)sp_len;
                int_fmin = int_sp;
                C_fmin = it->second;
            }
        }
        if (found){
            curr_substr = {len_fmin, int_fmin, C_fmin, start}; // {len_fmin, int_fmin, C_fmin, start};
            // still don't know if this is the correct fmin
                
            if (w_fmin > curr_substr){ // Compare LEXICOGRAPHICALLY Only strings of the same length will be compared
                all_fmin.clear();
                w_fmin = curr_substr;
            } else {
                while (all_fmin.back() > curr_substr) {
                    all_fmin.pop_back();
                }
            }
            all_fmin.push_back(curr_substr);

        }
        
        // Check if we are in a kmer
        //if (end - kmer_start + 1 == k){
        if (start > k-2){ // Once start reaches k-1 it means that we have looked for all the finimizers in the first k-mer
            
            while (get<3>(w_fmin) < kmer_start) {// {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                all_fmin.pop_front();
                w_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, uint64_t, int64_t, int64_t>{k+1,0,0,kmer_start};
            }
            
            if (all_fmin.size()>0){
                // This avoids storing the same finimizer multiple times for distinct kmers
                if (last_pos != get<3>(w_fmin) ) {Fmin.push_back(get<2>(w_fmin));} // pointer to C
                last_pos = get<3>(w_fmin);
            }
            
            kmer_start++;
        }
    }
    return Fmin;
}
   
   // FOUND COLORS: vector of bits and flip found and store in a vector
   // the flip back

   // BITMAPS SETS INSTEAD OF COLORS

    //TODO: int for the number of colors, change if needed
    void pseudoalignemnt_stats(vector<uint64_t>& Fmin, const vector<set<int>>& C, unordered_map<int, uint64_t>& results){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size();
        size_t rm_fmin = 0; // finimizers not found in the index
        
        //TODO if we knew the number of colors, results could be a vector of size colors

        for(const auto& f : Fmin){
            std::set<int> colors = C[f];//old hashTable.at(pair.first); // TODO replace this with colors you store the index of the finimizer (x) and then you can get the color from colors[x]
            for(const int& c : colors){
                results[c]++;
            }
        }
        return;
    }

    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const vector<set<int>>& C, vector<pair<int, float>>& results, const float& t){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size(); // # total finimizers
        size_t rm_fmin = 0; // finimizers not found in the index

        // count the number of colors found
        set<int> found_colors = {};
        std::unordered_map<int,uint64_t> fmin_per_color;
        
        for(const auto& f : Fmin){
            std::set<int> colors = C[f];
            for(const int& c : colors){ 
                found_colors.insert(c);
                fmin_per_color[c] ++;
            }
        }

        results.reserve(found_colors.size());

        for (const int& c : found_colors){         
            // For every color found, (#finimizers with that color)/(#tot finimizers - finimizers not found)
            float fraction = static_cast<float>(fmin_per_color[c]/static_cast<float>(found_fmin));
        
            if (t==1){ 
                if (fraction == t ){results.push_back({c,fraction});}
            } else{
                if (fraction > t){results.push_back({c,fraction});}
            }
        }
        return;
    }

//TODO remove?
// used in build
string print_finimizer_stats(const set<tuple<int64_t, int64_t, int64_t>>& finimizers, int64_t n_kmers, int64_t n_nodes, int64_t t){
    int64_t new_number_of_fmin = finimizers.size();
    int64_t sum_freq = 0;
    int64_t sum_len = 0;
    for (auto x : finimizers){
        sum_freq += get<1>(x);
        sum_len += get<0>(x);
    }

    string result = to_string(new_number_of_fmin) + "," + to_string(sum_freq) + "," + to_string(static_cast<float>(sum_freq) / static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len) / static_cast<float>(new_number_of_fmin)) + "," + to_string(n_kmers);

    write_log(to_string(t) + "," + result, LogLevel::MAJOR);
    write_log("#SBWT nodes: " + to_string(n_nodes) , LogLevel::MAJOR);
    write_log("#Distinct finimizers: " + to_string(new_number_of_fmin) , LogLevel::MAJOR);
    write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
    write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    return result;
}