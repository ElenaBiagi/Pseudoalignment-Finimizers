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
#include "BoundedDeque.hh"

#include "common.hh"
#include "bitsearch.hh"


// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? YES
// set ?
// TODO replace sbwt and LCS with const std::unordered_map<uint32_t, uint32_t>& B, const std::unordered_map<uint32_t, std::set<int>>& sB, const sdsl::bit_vector& T,
unordered_map<int64_t, uint64_t> rarest_fmin_streaming_search(const string& input, const unordered_map<uint32_t, pair<int64_t,int64_t> >& B, const std::unordered_map<uint32_t, int64_t>& sB, const std::vector<int_vector<1>>& T, const uint8_t plen, const int k){ 
    
    const int64_t str_len = input.size();

    unordered_map<int64_t, uint64_t> Fmin; // pointer to C, number of such finimizers
    Fmin.reserve(str_len-k+1);

    int64_t last_pos = 0;
    int64_t start = 0;
    int64_t kmer_start = 0; // start of the first k-mer

    bool found = false;
    uint8_t len_fmin = 0;
    //int64_t int_fmin = 0;
    string str_fmin;
    int64_t C_fmin = 0; // Can the result of color index be negative??? If not found??

    BoundedDeque<tuple<uint8_t, string, int64_t, int64_t>> all_fmin(input.size()-k+1);

    tuple<uint8_t, string, int64_t, int64_t> curr_substr; // length, fmin, C_offset, start
    tuple<uint8_t, string, int64_t, int64_t> w_fmin = {k+1,"1",0, input.size()}; // start will always be < str_len
    
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
        //const string prefix = input.substr(start,end);
        // Look for the prefix in B
        uint64_t int_p = prefix2int(input, start, plen);
        //if (B.find(intp) != B.end()) { // B contains all the possible prefixes of length p
        
        //int64_t pointer = B.find((uint32_t)int_p)->second;
        auto pp = B.at(int_p);
        int64_t pointer = pp.first;
        int64_t tails_so_far = pp.second;


        if (pointer != -1){
            // 2. Prefix found!
            string s = input.substr(start+plen, k-plen); // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch_new depending on tlen
            auto result = bitMagicSearch_new(T[pointer], s); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            //cerr << endl;
            
            if (result.first != -1){ 
                // b. Tail Found!
                found = true;
                len_fmin = plen + result.second; // TODO add the correct fmin len
                C_fmin = result.first + tails_so_far; // TODO NO NEED TO STORE THE COLORS NOW AS LONG AS WE KEEP THE OFFSET 
                //int_fmin = 0; // TODO add suffix () after prefix (int_p) // TODO strong the fmin as a string and not as a number 
                str_fmin = input.substr(start+plen, start+plen+result.second);
                reverse(str_fmin.begin(), str_fmin.end());
                // actual finimizer (as int) to check the colexicographically smallest one //TODO reverse it (NO, WRONG C,G) check 2 bits at a time?? but how to check simply? 
            }
        } else{
            // 1. Prefix NOT found
            // Shorten the prefix until you find a match
            // TODO should we keep the length of the SHORTEST finimizer? To know when to stop
            uint8_t sp_len = plen-1;
            uint64_t int_sp = prefix2int(input, start, sp_len); // It might be smarter to start from the longest prefix, done
            auto it = sB.find(int_sp);

            while(it == sB.end() and sp_len > 0){
                sp_len--;
                int_sp = prefix2int(input, start, sp_len);  // TODO add or remove one uint8_t at a time
            }
            if (sp_len != 0){
                found = true;
                len_fmin = (uint8_t)sp_len;
                //int_fmin = int_sp;
                str_fmin = input.substr(start,sp_len);
                // TODO this could be turned into numbers again now
                reverse(str_fmin.begin(), str_fmin.end());
                C_fmin = it->second;
            }
        }
        if (found){
            curr_substr = {len_fmin, str_fmin, C_fmin, start}; // {len_fmin, int_fmin, C_fmin, start};
            // still don't know if this is the correct fmin
                
            if (w_fmin > curr_substr){ // TODO compare colex easily // compare the REVERSE of strings and not int. Only strings of the same length will be compared
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
                w_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<uint8_t, string, int64_t, int64_t>{k+1,"0",0,kmer_start};
            }
            
            if (all_fmin.size()>0){
                // NO, we want to store them all -> This avoids storing the same finimizer (same) multiple times for distinct kmers
                //if (last_pos != get<3>(w_fmin) ) {Fmin.push_back(input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin)))};
                //last_pos = get<3>(w_fmin);
                // TODO improve this: e.g. store a counter for every finimizer
                Fmin[get<2>(w_fmin)]++;
                // sum of found finimizers = found k-mers
                //count++; // counts the number of kmers?? not used now
            }
            
            kmer_start++;
        }
    }
    return Fmin;
}
    

    /*
    unordered_map<string, uint64_t> rarest_fmin_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input){ 
    const int64_t n_nodes = sbwt.number_of_subsets();
    const int64_t k = sbwt.get_k();
    const vector<int64_t>& C = sbwt.get_C_array();
    BoundedDeque<tuple<int64_t, int64_t, int64_t, int64_t>> all_fmin(input.size());
    const int64_t str_len = input.size();
    tuple<int64_t, int64_t, int64_t, int64_t> w_fmin = {n_nodes,k+1,n_nodes,str_len+1}; // {freq, len, I start, end}

    int64_t freq;
    int64_t count = 0;
    int64_t start = 0;
    int64_t end;
    int64_t kmer_start = 0;
    pair<int64_t, int64_t> I = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new;
    int64_t I_start;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;

    unordered_map<string, uint64_t> Fmin;

    Fmin.reserve(str_len-k+1);
    int64_t last_pos = 0;
    
    // the idea is to start from the first pos which is i and move until finding something of ok freq
    // then drop the first char keeping track of which char you are starting from
    // Start is always < k as start <= end and end <k
    // if start == end than the frequency higher than t
    for (end = 0; end < str_len; end++) {
        char c = static_cast<char>(input[end] & ~32); // convert to uppercase using a bitwise operation //char c = toupper(input[i]);
        int64_t char_idx = get_char_idx(c);
        if (char_idx == -1) [[unlikely]]{
            cerr << "Error: unknown character: " << c << endl;
            cerr << "This works with the DNA alphabet = {A,C,G,T}" << endl;
            return {};
        } else {
            // 1) fmin interval
            I_new = sbwt.update_sbwt_interval(&c, 1, I);
            // (1) Finimizer(subseq) NOT found
            while(I_new.first == -1){
                //kmer_start = ++start;
                start++;
                if (start>end)[[unlikely]]{
                    I_new = {0, n_nodes - 1};
                    break;
                }
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = sbwt.update_sbwt_interval(&c, 1, I);
            }
            I = I_new;
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (2) Finimizer(subseq) freq > 0
            
            // (2b) Finimizer found
            if (freq ==1){ // 1. rarest
                while (freq == 1) { // 2. shortest
                    I_start = I.first;
                    curr_substr = {freq, end - start + 1, I_start, end};
                    // 2. drop the first char
                    // When you drop the first char you are sure to find x_2..m since you found x_1..m before
                    start ++;
                    I = drop_first_char(end - start + 1, I, LCS, n_nodes);
                    freq = (I.second - I.first + 1);
                }
                if (w_fmin > curr_substr) {
                    all_fmin.clear();
                    w_fmin = curr_substr;
                } else{
                    while (all_fmin.back() > curr_substr) {
                        all_fmin.pop_back();
                    }
                }
                all_fmin.push_back(curr_substr);
            }

            // Check if we are in a kmer
            if (end - kmer_start + 1 == k){
            
                count++; // counts the number of kmers?? not used now
                while (((get<3>(w_fmin)-get<1>(w_fmin)) +1) < kmer_start) {// {freq, len, I start, end}
                    all_fmin.pop_front();
                    w_fmin = (all_fmin.size()>0) ? all_fmin.front() : tuple<int64_t, int64_t, int64_t, int64_t>{n_nodes,k+1,n_nodes,kmer_start+k};
                }
                
                if (all_fmin.size()>0){
                    // NO, we want to store them all -> This avoids storing the same finimizer (same) multiple times for distinct kmers
                    //if (last_pos != get<3>(w_fmin) ) {Fmin.push_back(input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin)))};
                    //last_pos = get<3>(w_fmin);
                    // TODO improve this: e.g. store a counter for every finimizer
                    string finimizer = input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin));
                    Fmin[finimizer]++;
                }
                

                kmer_start++;

                //I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
            }
        }
    }
    //if (count != Fmin.size()){
    std::cerr << "total k-mers = " << count << ", total finimizers = " << Fmin.size()<< std::endl;
    //}
    return Fmin;
}



    void pseudoalignemnt_stats(unordered_map<string, uint64_t>& Fmin, const std::unordered_map<std::string,std::set<int>>& hashTable, unordered_map<int, uint64_t>& results){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size();
        size_t rm_fmin = 0; // finimizers not found in the index

        // count the number of colors found
        set<int> found_colors = {};
        vector<set<int>> found_colors_single = {};
        found_colors_single.reserve(found_fmin);
        
        for(const auto& pair : Fmin){
            try {
                std::set<int> colors = hashTable.at(pair.first);
                found_colors_single.push_back(colors);
                for(const int& c : colors){
                    found_colors.insert(c);
                    results[c]+=pair.second;
                }
            } catch (const std::out_of_range& e) {
                rm_fmin+=pair.second; // this is not used
            }
        }
        //if ( rm_fmin > 0 ) std::cerr << rm_fmin << std::endl;
        
        //std::cerr << found_colors.size() << " found colors" << std::endl;
        return;
    }

     void pseudoalignemnt_stats(unordered_map<string, uint64_t>& Fmin, const std::unordered_map<std::string,std::set<int>>& hashTable, vector<pair<int, float>>& results, const float& t){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size(); // # distinct finimizers
        size_t ok_fmin = 0; // total finimizers found
        size_t rm_fmin = 0; // finimizers not found in the index

        // count the number of colors found
        set<int> found_colors = {};
        std::unordered_map<int,uint64_t> fmin_per_color;
        
        for(const auto& pair : Fmin){
            try {
                std::set<int> colors = hashTable.at(pair.first);
                ok_fmin += pair.second;
                for(const int& c : colors){ 
                    found_colors.insert(c);
                    fmin_per_color[c] += pair.second;
                 }
            } catch (const std::out_of_range& e) {
                rm_fmin+=pair.second; // not used now
            }
        }
        //if ( rm_fmin > 0 ) std::cerr << rm_fmin << std::endl;
        
        // TODO all colors and not only the found ones (input?)
        results.reserve(found_colors.size());

        for (const int& c : found_colors){
            
            //std::cerr<< "color " << c << ": " << fmin_per_color[c] << " finimizers" << std::endl;
           
            // For every color found, (#finimizers with that color)/(#tot finimizers - finimizers not found)
            float fraction = static_cast<float>(fmin_per_color[c]/static_cast<float>(ok_fmin));
        
            if (t==1){ 
                if (fraction == t ){results.push_back({c,fraction});}
            } else{
                if (fraction > t){results.push_back({c,fraction});}
            }
        
        }
        return;
    } */
   
    void pseudoalignemnt_stats(unordered_map<int64_t, uint64_t>& Fmin, const vector<set<int>>& C, unordered_map<int, uint64_t>& results){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size();
        size_t rm_fmin = 0; // finimizers not found in the index

        // count the number of colors found
        set<int> found_colors = {};
        vector<set<int>> found_colors_single = {};
        found_colors_single.reserve(found_fmin);
        
        for(const auto& pair : Fmin){
            try {
                std::set<int> colors = C[pair.first];//old hashTable.at(pair.first); // TODO replace this with colors you store the index of the finimizer (x) and then you can get the color from colors[x]
                found_colors_single.push_back(colors);
                for(const int& c : colors){
                    found_colors.insert(c);
                    results[c]+=pair.second;
                }
            } catch (const std::out_of_range& e) {
                rm_fmin+=pair.second; // this is not used
            }
        }
        //if ( rm_fmin > 0 ) std::cerr << rm_fmin << std::endl;
        
        //std::cerr << found_colors.size() << " found colors" << std::endl;
        return;
    }

    void pseudoalignemnt_stats(unordered_map<int64_t, uint64_t>& Fmin, const vector<set<int>>& C, vector<pair<int, float>>& results, const float& t){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size(); // # distinct finimizers
        size_t ok_fmin = 0; // total finimizers found
        size_t rm_fmin = 0; // finimizers not found in the index

        // count the number of colors found
        set<int> found_colors = {};
        std::unordered_map<int,uint64_t> fmin_per_color;
        
        for(const auto& pair : Fmin){
            try {
                std::set<int> colors = C[pair.first];
                ok_fmin += pair.second;
                for(const int& c : colors){ 
                    found_colors.insert(c);
                    fmin_per_color[c] += pair.second;
                 }
            } catch (const std::out_of_range& e) {
                rm_fmin+=pair.second; // not used now
            }
        }
        //if ( rm_fmin > 0 ) std::cerr << rm_fmin << std::endl;
        
        // TODO all colors and not only the found ones (input?)
        results.reserve(found_colors.size());

        for (const int& c : found_colors){
            
            //std::cerr<< "color " << c << ": " << fmin_per_color[c] << " finimizers" << std::endl;
           
            // For every color found, (#finimizers with that color)/(#tot finimizers - finimizers not found)
            float fraction = static_cast<float>(fmin_per_color[c]/static_cast<float>(ok_fmin));
        
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
 // ????
/*     std::ostream& operator<<(std::ostream& os, const std::set<int>& set) {
    os << "{";
    for (auto it = set.begin(); it != set.end(); ++it) {
        os << *it;
        if (std::next(it) != set.end()) {
            os << ", ";
        }
    }
    os << "}";
    return os;
} */

