#pragma once

#include <string>
#include <cstring>
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


// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? NO
// set ?
vector<string> rarest_fmin_streaming_search(const plain_matrix_sbwt_t& sbwt, const sdsl::int_vector<>& LCS, const string& input){ 
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
    pair<int64_t, int64_t> I = {0, n_nodes - 1}, I_kmer = {0, n_nodes - 1};
    pair<int64_t, int64_t> I_new, I_kmer_new;
    int64_t I_start;
    tuple<int64_t, int64_t, int64_t, int64_t> curr_substr;

    vector<string> Fmin;
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
                kmer_start = ++start;
                if (start>end)[[unlikely]]{
                    I_new = {0, n_nodes - 1};
                    I_kmer = I_new;
                    break;
                }
                I = drop_first_char(end - start, I, LCS, n_nodes); // The result (substr(start++,end)) cannot have freq == 1 as substring(start,end) has freq >1
                I_new = sbwt.update_sbwt_interval(&c, 1, I);
                I_kmer = I_new;
            }
            I = I_new;
            freq = (I.second - I.first + 1);
            I_start = I.first;
            // (2) Finimizer(subseq) freq > 0
            // Check if the Kmer interval has to be updated
            if ( start != kmer_start){
                I_kmer_new = sbwt.update_sbwt_interval(&c, 1, I_kmer);
                while(I_kmer_new.first == -1){
                    // kmer NOT found
                    kmer_start++;
                    I_kmer = drop_first_char(end - kmer_start, I_kmer, LCS, n_nodes);
                    I_kmer_new = sbwt.update_sbwt_interval(&c, 1, I_kmer);
                } 
                I_kmer = I_kmer_new;
            } else { 
                I_kmer = I;
            }
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
            
            //TODO do we want to keep this check??
            // Check if the kmer is found
            if (end - kmer_start + 1 == k){
            
                count++; // counts the number of kmers?? not used now
                while ((get<3>(w_fmin)-get<1>(w_fmin) +1) < kmer_start) {
                    all_fmin.pop_front();
                    w_fmin = all_fmin.front();
                }
                
                if (last_pos != get<3>(w_fmin) ) Fmin.push_back(input.substr(get<3>(w_fmin)-get<1>(w_fmin)+1,get<1>(w_fmin)));
                last_pos = get<3>(w_fmin);

                kmer_start++;
                I_kmer = drop_first_char(end - kmer_start + 1, I_kmer, LCS, n_nodes);
            }
        }
    }
    return Fmin;
}

    void pseudoalignemnt_stats(vector<string>& Fmin, const std::unordered_map<std::string,std::set<int>>& hashTable, vector<pair<int,float>>& results, set<int>& intersection){
        // count the number of finimizers found
        size_t found_fmin = Fmin.size();
        
        // count the number of colors found
        set<int> found_colors = {};
        vector<set<int>> found_colors_single = {};
        found_colors_single.reserve(found_fmin);
        for(const string& fmin : Fmin){
            //std::cerr << fmin << ", ";
            set<int> colors = hashTable.at(fmin);
            found_colors_single.push_back(colors);
            for(const int& c : colors){found_colors.insert(c); }
        }
        //std::cerr << std::endl;
        //std::cerr << found_fmin << " found Finimizers" << std::endl;
        //std::cerr << found_colors.size() << " found colors" << std::endl;

        // for every color found, (#finimizers with that color)/(#tot finimizers)
        for (const int& c : found_colors){
            int64_t c_found_fmin = 0;
            for (const set<int>& f : found_colors_single){
                if (f.count(c)){ c_found_fmin++;}
            }
            float fraction = static_cast<float>(c_found_fmin/static_cast<float>(found_fmin));
            if (fraction > 0.8){results.push_back({c,fraction});}
            //std::cerr<< "color " << c << ": " << c_found_fmin << " finimizers" << std::endl;
            //std::cerr << c << "; " << static_cast<float>(c_found_fmin/static_cast<float>(found_fmin)) << std::endl;
        }
            //for (const std::pair<int, float>& p : results) { std::cerr << "{ " << p.first << ", " << p.second << " } " << std::endl; }
        // intersection
        intersection = found_colors_single[0];

        // Iterate through the remaining sets
        for (size_t i = 1; i < found_colors_single.size(); ++i) {
            set<int> temp;
            set_intersection(intersection.begin(), intersection.end(),
                                found_colors_single[i].begin(), found_colors_single[i].end(),
                                std::inserter(temp, temp.begin()));
            intersection = std::move(temp);
        }
        
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

    while ((end = input.find('N', start)) != std::string::npos) {
        if (end > start+k-2) { // Exclude strings shorter than k
            result.emplace_back(input.substr(start, end - start));
        }
        start = end + 1;
    }

    // Add the last part if it's at least k characters long
    if (start+k-2< input.size()) {
        result.emplace_back(input.substr(start));
    }

    return result;
}