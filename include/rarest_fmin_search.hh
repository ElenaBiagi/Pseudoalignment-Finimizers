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
#include "Buckets.hh"


// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? YES
// set ?
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 
    const int64_t str_len = input.size();

    vector<uint64_t> Fmin;// pointer to C

    //int64_t last_pos = 0;
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
        
        //int64_t tails_so_far = B[int_p];

        if (buckets[int_p].has_value()){
            // 2. Prefix found!
            const Bucket& bucket_p = *buckets[int_p];
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch_new depending on tlen

            char s_len = (str_len >= start+k) ? k-plen : str_len-start;
            s_int = prefix2int(input, start+plen, s_len); // tail
            auto result = bitMagicSearch_new(bucket_p.tail_data, s_int, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            
            if (result.first != -1){ 
                // b. Tail Found!
                found = true;
                len_fmin = plen + result.second; // TODO add the correct fmin len
                C_fmin = bucket_p.color_set_ids[result.first]; // TODO NO NEED TO STORE THE COLORS NOW AS LONG AS WE KEEP THE OFFSET 
                // Store the string as a number. OK as only strings of the same length will be compared 
                
                // Extract the number instead of converting again
                int_fmin = (s_int >> ((s_len - len_fmin) * 2)) & ((1ULL << (len_fmin * 2)) - 1); 
                //old //int_fmin =  prefix2int(input, start+plen, len_fmin);                
            }
        } else{
            // 1. Prefix NOT found
            // Start from the shortest possible prefix
            // if you find a real match, stop
            uint8_t sp_len = 1;
            uint64_t int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1); // prefix of length 1 of int_p
            auto it = sB.find(int_sp);

            while (sp_len < plen){ // sp_len must be < plen as the whole prefix was not found
                if (it != sB.end() && sp_len == it->second.first){ // real match
                    found = true;
                    len_fmin = sp_len;
                    int_fmin = int_sp;
                    C_fmin = it->second.second;
                    break;
                }
                sp_len++;
                int_sp = (int_p >> ((plen - sp_len) * 2)) & ((1ULL << (sp_len * 2)) - 1);
                it = sB.find(int_sp);
            }
        }
        if (found){
            curr_substr = {len_fmin, int_fmin, C_fmin, start}; // {len_fmin, int_fmin, C_fmin, start};
            // still don't know if this is the correct fmin
                
            if (w_fmin > curr_substr){ // Compare LEXICOGRAPHICALLY Only strings of the same length will be compared so it's ok to use numbers
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
            
            if (all_fmin.size()>0){Fmin.push_back(get<2>(w_fmin));}
            
            kmer_start++;
        }
    }
    return Fmin;
}
   
   // TODO: FOUND COLORS: vector of bits and flip found and store in a vector, then flip back
   // BITMAPS SETS INSTEAD OF COLORS

    //TODO: int for the number of colors, change if needed
    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<uint64_t>& results){ 
        // count the number of finimizers found
        results.assign(n_colors, 0);
        
        #pragma omp parallel for
        for(const auto& start : Fmin){
            for (uint64_t i=0; i< n_colors; i++){
                #pragma omp atomic
                results[i]+=color_sets_concat[start+i];
            }
        }
        return;
    }

    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size(); // # total finimizers
        //size_t rm_fmin = 0; // finimizers not found in the index

        
        vector<uint64_t> tot_res;
        tot_res.assign(n_colors, 0);

        #pragma omp parallel for
        for(const auto& start : Fmin){
            for (uint64_t i=0; i< n_colors; i++){
                #pragma omp atomic
                tot_res[i]+=color_sets_concat[start+i];
            }
        }
        results.assign(n_colors, 0);
        
        #pragma omp parallel for
        for (uint64_t i=0; i< n_colors; i++){         
            // For every color found, (#finimizers with that color)/(#tot finimizers)
            float fraction = static_cast<float>(tot_res[i]/static_cast<float>(found_fmin));
            
            if (t==1 && fraction ==1){ 
                results[i]=1;
            } else if (fraction > t){
                results[i]=fraction;
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