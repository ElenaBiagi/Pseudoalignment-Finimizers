#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>


//#include "SeqIO.hh"
#include "BoundedDeque.hh"
#include "CircularBuffer.hh"

#include "common.hh"
#include "bitsearch.hh"
#include "Buckets.hh"


void FindShortFinimizer(const uint64_t int_sp_len, uint64_t int_sp, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint64_t start,const uint64_t end, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t>& k_fmin){
    // 2. Prefix NOT found
    // Start from the longest possible prefix
    // if you find a real match, stop
    uint64_t sp_len = int_sp_len; // plen is here plen-1: sp_len must be < plen as the whole prefix was not found  
    while (sp_len > 0){ 
        auto it = sB.find(int_sp);
        if (it != sB.end() && sp_len == it->second.first){ // real match 
            //all_fmin.insert(make_tuple(sp_len, int_sp, it->second.second, start));
            if ((start + sp_len -1) > end) { next_candidates.push_back(make_tuple(sp_len+start-1, int_sp, it->second.second, start));} // Sorted based on END
            else { 
                tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {sp_len, int_sp, it->second.second, start};
                //tuple<uint64_t, uint64_t, uint64_t, uint64_t> k_fmin = curr_candidates.front();
                if (new_fmin < k_fmin) {
                    curr_candidates.clear();
                    k_fmin = new_fmin;
                }
                else{
                    while (curr_candidates.back()> new_fmin) {curr_candidates.pop_back();}
                }
                curr_candidates.push_back(new_fmin); 
            }
            return;
        }
        int_sp >>= 2;
        sp_len--;
    }
}

void FindPrefix(const vector<optional<Bucket>>& buckets, const uint64_t plen, const char s_len, const uint64_t int_s, const uint64_t int_p, const uint64_t start, const uint64_t end, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t>& k_fmin){
    const Bucket& bucket_p = *buckets[int_p];
    auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len);  
    if (pos > -1){
        uint64_t f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); //  Shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        if (plen+len ==32){cerr<< "ERROR plen+len " << plen << " + "<< len<< endl; }
        if ((start + plen+len -1) > end) { next_candidates.push_back(make_tuple(plen+len+start-1, f_int, bucket_p.color_set_ids[pos], start));} // Sorted based on END
        else { 
            tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {plen + len, f_int, bucket_p.color_set_ids[pos], start};
            if (new_fmin < k_fmin) {
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else{
                while (curr_candidates.back()> new_fmin) {curr_candidates.pop_back();}
            }
            curr_candidates.push_back(new_fmin); 
        }
    }
}

// The query is shorter than (k - plen)
void FindPrefix_short(const vector<optional<Bucket>>& buckets, const uint64_t plen, const char s_len, const uint64_t int_s, const uint64_t int_p, const uint64_t start, const uint64_t end, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t>& k_fmin){
    const Bucket& bucket_p = *buckets[int_p];
    auto [pos,len] = bitMagicSearch_short(bucket_p.tail_data, int_s, s_len);  
    if (pos > -1){
        if (plen+len ==32){cerr<< "ERROR short plen+len " << plen << " + "<< len<< endl; }
        uint64_t f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); //  Shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        if ((start + plen+len -1) > end) { next_candidates.push_back(make_tuple(plen+len+start-1, f_int, bucket_p.color_set_ids[pos], start));} // Sorted based on END
        else { 
            tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {plen + len, f_int, bucket_p.color_set_ids[pos], start};
            if (new_fmin < k_fmin) {
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else{
                while (curr_candidates.back()> new_fmin) {curr_candidates.pop_back();}
            }
            curr_candidates.push_back(new_fmin); 
        }
    }
}

void PickFinimizer(vector<uint64_t>& Fmin, const uint64_t kmer_start, const uint64_t k, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t>& k_fmin){//, const string& input){
    //size_t str_len = input.size();
    if (!curr_candidates.empty()){
        k_fmin= curr_candidates.front();
        /* if (!Fmin.empty()){
            if( get<2>(k_fmin) == Fmin[Fmin.size()-1]){
            cout;}
        } else {
        cout << input.substr(get<3>(k_fmin),get<0>(k_fmin)) << endl;
        if (get<0>(k_fmin)==32){
            cerr << "ERROR LENGTH =32!!!" << endl;
        }
        }  */
        Fmin.push_back(get<2>(k_fmin));
    } else{
        k_fmin = static_cast<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>(make_tuple(k+1,0,0,kmer_start+1));
        cerr << "finimizer not found for kmer " << kmer_start << "-"<< kmer_start + k-1  << endl;// " " << input.substr(kmer_start, min((uint64_t)k, str_len - kmer_start)) << endl;
    }

    // 1. Check if this finimizer is good for the next k-mer (still in the window)
    while (!curr_candidates.empty() && get<3>(curr_candidates.front()) <= kmer_start){
        curr_candidates.pop_front();
    }
    k_fmin = (curr_candidates.empty()) ? static_cast<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>(make_tuple(k+1,0,0,kmer_start+1)) : curr_candidates.front();

    // 2. Check if the NEXT finimizer would be good for the next k-mer
    if (!next_candidates.empty()){
        const auto& next_fmin = next_candidates.front(); // tuple<uint64_t, uint64_t, uint64_t, uint64_t>
        tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {get<0>(next_fmin)-get<3>(next_fmin)+1, get<1>(next_fmin), get<2>(next_fmin), get<3>(next_fmin)};
        if (get<0>(next_fmin) <= kmer_start + k){ // end of the fmin is before end of next kmer
            if (new_fmin < k_fmin) { // always true if curr_candidates is empty
                curr_candidates.clear();
                k_fmin = new_fmin;
            }
            else{
                while (curr_candidates.back()> new_fmin) {curr_candidates.pop_back();}
            }
            curr_candidates.push_back(new_fmin);
            next_candidates.pop_front();
        }
    }
}
void rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint64_t plen, const uint64_t k, vector<uint64_t>& Fmin){ 

    const int64_t str_len = input.size();

    uint64_t start = 0;
    uint64_t kmer_start = 0;

    BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> curr_candidates(k); // sort based on len, int (color, start)
    BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>> next_candidates(k); // sort by end (start+len-1)
    tuple<uint64_t, uint64_t, uint64_t, uint64_t> k_fmin = {k+1,0,0,kmer_start};
    curr_candidates.push_back(k_fmin);
    
    uint64_t int_p = prefix2int(input, start, plen); // start = 0
    char s_len = k-plen;
    uint64_t int_s = prefix2int(input, start+plen, s_len); // tail
    uint64_t int_sp;

    // Check the first characters
    if (buckets[int_p].has_value()){
        // 1. prefix found
        FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input); 
    }else{
        // 2. look for a shorter finimizer
        int_sp = int_p >> 2;
        FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input); 
    }
    
    // The first k-1 characters do not contail all possible finimizers for the first k-mer
    for (start = 1; start < k-1; start++ ){ 
        int_p = stream_kmer(int_p, input[start + plen - 1], plen); // shorten by 1 at every loop iteration 
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
        int_s = prefix2int(input, start+plen, s_len);; // tail   
        if (buckets[int_p].has_value()){
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input); 
        }
        else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input); 
        }
    }

    uint64_t ss = start;
    s_len = k - plen; // Constant s_len
    for (start = ss ; start <= str_len-k; start++ ){ 
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen 
        int_s = stream_kmer(int_s, input[start + plen + s_len - 1], s_len);
        
        if (buckets[int_p].has_value()){
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input);  
        }else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input);  
        }
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); // , input); 
        kmer_start++;
    }
    // Shorter s_len
    ss = start;
    for (start = ss ; start < str_len-plen+1; start++ ){ 
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        s_len = str_len-start-plen;
        int_s &= ((1ULL << (2 * s_len)) - 1); // Shorten int_s by 2 at the beginning
        
        if (buckets[int_p].has_value()){        
            FindPrefix_short(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input);  
        }else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input);  
        }
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); // , input); 
        kmer_start++;
    }

    //The last plen-1 characters cannot contain a prefix
    ss = start;
    uint64_t s_plen = plen;
    for (start = ss; start < str_len; start++ ){ 
        s_plen--;
        int_p &= ((1ULL << (2 * s_plen)) - 1); // Shorten int_p by 2 
        FindShortFinimizer(s_plen, int_p, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // , input); // this shortens s_plen by 1 internally
        
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); // , input); 
        kmer_start++;
    }
    return;
}

void read_colors(const uint64_t* data, const uint64_t n_colors, vector<uint64_t>& tot_res, const vector<pair<uint64_t, uint64_t>> fmin_v){
    for (const auto& [start,freq] : fmin_v) {
        const uint64_t* ptr = data + (start * n_colors) / 64;
        uint64_t bit_offset = (start * n_colors) % 64;

        uint64_t color_id = 0;

        // 1. Read the first word
        uint64_t bits_to_read = std::min(64UL - bit_offset, n_colors);
        uint64_t mask = (bits_to_read == 64) ? ~0ULL : ((1ULL << bits_to_read) - 1);
        uint64_t word = (*ptr >> bit_offset) & mask;

        while (word != 0) {
            uint64_t bit = __builtin_ctzll(word);
            tot_res[bit]+=freq;
            word &= word - 1;
        }

        ++ptr;
        color_id += bits_to_read;

        // 2. Read aligned words in btw
        while (color_id + 64 <= n_colors) {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;) {
                uint64_t bit = __builtin_ctzll(w);
                tot_res[color_id + bit]+=freq;
                w &= w - 1;
            }
            color_id += 64;
        }

        // 3. Read the last word (if any)
        uint64_t bits_left = n_colors - color_id;
        if (bits_left > 0) {
            uint64_t mask = ((1ULL << bits_left) - 1);
            uint64_t word = *ptr & mask;

            while (word != 0) {
                uint64_t bit = __builtin_ctzll(word);
                tot_res[color_id + bit]+=freq;
                word &= word - 1;
            }
        }
    }
}

//slower than counting sort
void bounded_counting_sort (const std::vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors, const uint64_t min_value){    
    // the range of ok values goes from min_value to found_fmin

    vector<uint16_t> counts;
    counts.resize(found_fmin+1);

    for (size_t idx = 0; idx < n_colors; idx++) {
        counts[results[idx]]++;
    }

    // Cumulative Sums
    for (size_t c = 1; c < counts.size(); c++) {
        counts[c]+=counts[c-1];
    }

    ans.resize(n_colors);
    for (size_t idx = 0; idx < n_colors; idx++) {
        const int64_t new_idx = results[idx] - min_value;
        if (new_idx >=0) {
            ans[counts[results[idx]] - 1] = {static_cast<uint16_t>(idx), results[idx]};
            counts[results[idx]]--;
        } 
    }
}

//slower than counting sort
void new_bounded_counting_sort (const std::vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors, const uint64_t min_value){    
    // Store only the colors with #matches >= min_value

    vector<uint16_t> counts;
    counts.resize(found_fmin-min_value+1);

    for (size_t idx = 0; idx < n_colors; idx++) {
        if (results[idx] >= min_value) {
            const int64_t new_idx = results[idx] - min_value;
            counts[new_idx]++;
        }
    }

    // Cumulative Sums
    for (size_t c = 1; c < counts.size(); c++) {
        counts[c]+=counts[c-1];
    }

    ans.resize(counts.back()); // only store some colors
    for (size_t idx = 0; idx < n_colors; idx++) {

        const int64_t new_idx = results[idx] - min_value;
        if (results[idx]>=min_value) {
            ans[counts[new_idx] - 1] = {static_cast<uint16_t>(idx), results[idx]};
            counts[new_idx]--;
        } 
    }
}

//slower than counting sort
void newnew_bounded_counting_sort (const std::vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors, const uint64_t min_value){    
    // Store only the colors with #matches >= min_value

    vector<uint16_t> counts;
    counts.resize(found_fmin+1);

    for (size_t idx = 0; idx < n_colors; idx++) {
        counts[results[idx]]++;
    }

    // Cumulative Sums
    for (size_t c = min_value +1; c < counts.size(); c++) {
        counts[c]+=counts[c-1];
    }

    ans.resize(counts.back()); // only store some colors
    for (size_t idx = 0; idx < n_colors; idx++) {
        if (results[idx]>=min_value) {
            ans[counts[results[idx]] - 1] = {static_cast<uint16_t>(idx), results[idx]};
            counts[results[idx]]--;
        } 
    }
}


void counting_sort (const std::vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors){    
    vector<uint16_t> counts;
    counts.resize(found_fmin+1); 
    for (size_t idx = 0; idx < n_colors; idx++) {
        counts[results[idx]]++;
    }

    // Cumulative Sums
    for (size_t c = 1; c < counts.size(); c++) {
        counts[c]+=counts[c-1];
    }

    ans.resize(n_colors);
    for (size_t idx = 0; idx < n_colors; idx++) {
        ans[counts[results[idx]] - 1] = {static_cast<uint16_t>(idx), results[idx]};
        counts[results[idx]]--;
    }
}


void pseudoalignment_stats(vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
    vector<uint64_t> results;
    results.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();

    // Count freq of each fmin
    std::unordered_map<uint64_t, uint64_t> fmin_counts;
    fmin_counts.reserve(Fmin.size());
    for (const auto& v: Fmin) {
        fmin_counts[v]++;
    }

    // vector for sorted output so that it is possible to scan color_set_concat ????
    std::vector<std::pair<uint64_t, uint64_t>> fmin_v(fmin_counts.begin(), fmin_counts.end());
    std::sort(fmin_v.begin(), fmin_v.end()); 
    read_colors(data, n_colors, results, fmin_v);

    const size_t found_fmin = Fmin.size(); // # total finimizers
    counting_sort(results, ans, found_fmin, n_colors);
    return;
}

uint16_t pseudoalignment_stats(vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float& t){ 
    vector<uint64_t> tot_res;
    tot_res.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();
    
    // Count freq of each fmin
    std::unordered_map<uint64_t, uint64_t> fmin_counts;
    for (auto v : Fmin) {
        fmin_counts[v]++;
    }

    // vector for sorted output so that it is possible to scan color_set_concat

    std::vector<std::pair<uint64_t, uint64_t>> fmin_v(fmin_counts.begin(), fmin_counts.end());
    std::sort(fmin_v.begin(), fmin_v.end());

    read_colors(data, n_colors, tot_res, fmin_v);

    // Check the values above the minimum
    const size_t found_fmin = Fmin.size(); // # total finimizers
    const uint64_t min_value = found_fmin * t;

    // TODO IMPROVE
    //newnew_bounded_counting_sort(tot_res, ans, found_fmin, n_colors, min_value);
    counting_sort(tot_res, ans, found_fmin, n_colors);

    /* vector<uint64_t> results;
    results.reserve(n_colors);
    // TODO DO NOT REMOVE VALUES HERE YET
    for (uint64_t i = 0; i < n_colors; ++i) {
        if (tot_res[i] > min_value) {
            results.emplace_back(i, tot_res[i]);
        }
    } */

    
    return min_value;
}
