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

char get_char_idx(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

/* 
TODO REMOVE
// Returns the end point (inclusice) of the first k-mer in the concatenation of the unitigs
int64_t lookup_from_branch_dictionary(int64_t kmer_colex, int64_t k, const sdsl::rank_support_v5<>& Ustart_rs, const PackedStrings& unitigs){
    int64_t unitig_rank = Ustart_rs.rank(kmer_colex);
    assert(unitig_rank < unitigs.ends.size());
    int64_t global_unitig_start = 0;
    if(unitig_rank > 0) global_unitig_start = unitigs.ends[unitig_rank-1];
    return global_unitig_start + k - 1;
}

int64_t lookup_from_finimizer_dictionary(int64_t finimizer_colex, const sdsl::rank_support_v5<>& fmin_rs, const sdsl::int_vector<>& global_offsets){
    int64_t finimizer_id = fmin_rs.rank(finimizer_colex);
    return global_offsets[finimizer_id];
}
 */

// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? YES
// set ?
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

    void pseudoalignemnt_stats(unordered_map<string, uint64_t>& Fmin, const std::unordered_map<std::string,std::set<int>>& hashTable, unordered_map<int, uint64_t>& results, set<int>& intersection, const float& t){ // vector<pair<int,float>>& results
        // count the number of finimizers found
        size_t found_fmin = Fmin.size();
        size_t rm_fmin = 0; // finimizers not found in the index

        // count the number of colors found
        set<int> found_colors = {};
        vector<set<int>> found_colors_single = {};
        found_colors_single.reserve(found_fmin);
        
        for(const auto& pair : Fmin){
            //std::cerr << fmin << ", ";
            try {
                std::set<int> colors = hashTable.at(pair.first);
                found_colors_single.push_back(colors);
                for(const int& c : colors){
                    found_colors.insert(c);
                    results[c]+=pair.second;
                }
            } catch (const std::out_of_range& e) {
                rm_fmin+=pair.second;
            }
        }

        //if ( rm_fmin > 0 ) std::cerr << rm_fmin << std::endl;
        
        //std::cerr << found_fmin - rm_fmin << " found Finimizers" << std::endl;
        //std::cerr << found_colors.size() << " found colors" << std::endl;

        // TODO all colors and not only the found ones (input?)
        // it exists already vector<std::pair<int,uint64_t>> results;
        results.reserve(found_colors.size());

        /* for (const int& c : found_colors){
            int64_t c_found_fmin = 0;
            for (const set<int>& f : found_colors_single){
                if (f.count(c)){c_found_fmin++;}
            }
            std::pair<int,uint64_t>  p = {c,c_found_fmin};
            results.emplace_back(p);
            //std::cerr<< "color " << c << ": " << c_found_fmin << " finimizers" << std::endl;

 */        /*    
            // for every color found, (#finimizers with that color)/(#tot finimizers - finimizers not found)


            float fraction = static_cast<float>(c_found_fmin/static_cast<float>(found_fmin-rm_fmin));
        
        
            // TODO FIX INTERSECTION
            if (t==1){ 
                if (fraction == t ){results.push_back({c,fraction});}
            } else{
                if (fraction > t){results.push_back({c,fraction});}
            } */
        
            //std::cerr << c << "; " << static_cast<float>(c_found_fmin/static_cast<float>(found_fmin)) << std::endl;
        //}
            //for (const std::pair<int, float>& p : results) { std::cerr << "{ " << p.first << ", " << p.second << " } " << std::endl; }
        
        // intersection
        /* 
        if (found_colors_single.empty()){
            intersection = {};
        } else {
        
            intersection = found_colors_single[0];

            // Iterate through the remaining sets
            for (size_t i = 1; i < found_colors_single.size(); ++i) {
                set<int> temp;
                set_intersection(intersection.begin(), intersection.end(),
                                    found_colors_single[i].begin(), found_colors_single[i].end(),
                                    std::inserter(temp, temp.begin()));
                intersection = std::move(temp);
            } 
        }*/
        
        return;
    }


//TODO remove?
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

    void printHashTable(const std::unordered_map<std::string, std::set<int>>& hashTable) {
        std::cerr << "HASH TABLE" << std::endl;
        for (const auto& pair : hashTable) {
            std::cerr << "Key: " << pair.first << ", Values: " << pair.second << std::endl;
        }
        std::cerr << "HASH TABLE done" << std::endl;

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
    result.erase(std::remove(result.begin(), result.end(), 'N'), result.end());
    return result; // Return the resulting string as const
}

vector< std::string> split_by_N(const std::string &input, const int64_t k) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t end = 0;

    while ((end = input.find_first_of("BDEFHIJKLMNOPQRSUVWXYZbdefhijklmnopqrstuvwxyz", start)) != std::string::npos) {
        if (end - start +1 >= k) { // Exclude strings shorter than k
            result.emplace_back(input.substr(start, end - start));
        }
        start = end + 1;
    }

    // Add the last part if it's at least k characters long
    if (input.size()-start >=k) {
        result.emplace_back(input.substr(start));
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