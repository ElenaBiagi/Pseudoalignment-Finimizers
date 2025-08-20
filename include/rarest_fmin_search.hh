#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

#include "BoundedDeque.hh"

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

void PickFinimizer(vector<int64_t>& Fmin, const uint64_t kmer_start, const uint64_t k, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t>& k_fmin){//, const string& input){
    if (!curr_candidates.empty()){
        k_fmin= curr_candidates.front();
        
        // cout << input.substr(get<3>(k_fmin),get<0>(k_fmin)) << endl;
        Fmin.push_back(get<2>(k_fmin));
    } else{
        Fmin.push_back(-1); // This ensures that Fmin and r_Fmin have the same length
        //k_fmin = static_cast<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>(make_tuple(k+1,0,0,kmer_start+1));
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
void rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint64_t plen, const uint64_t k, vector<int64_t>& Fmin){ 

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

void read_colors(const uint64_t* data, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<int64_t, uint64_t>>& fmin_v){
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
            results[bit]+=freq;
            word &= word - 1;
        }

        ptr++;
        color_id += bits_to_read;

        // 2. Read aligned words in btw
        while (color_id + 64 <= n_colors) {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;) {
                uint64_t bit = __builtin_ctzll(w);
                results[color_id + bit]+=freq;
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
                results[color_id + bit]+=freq;
                word &= word - 1;
            }
        }
    }
}

void read_colors_old(const uint64_t* data, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<uint64_t, uint64_t>>& fmin_v){
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
            results[bit]+=freq;
            word &= word - 1;
        }

        ++ptr;
        color_id += bits_to_read;

        // 2. Read aligned words in btw
        while (color_id + 64 <= n_colors) {
            uint64_t word = *ptr++;
            for (uint64_t w = word; w != 0;) {
                uint64_t bit = __builtin_ctzll(w);
                results[color_id + bit]+=freq;
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
                results[color_id + bit]+=freq;
                word &= word - 1;
            }
        }
    }
}

//slower than counting sort
void bounded_counting_sort (const vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors, const uint64_t min_value){    
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
void new_bounded_counting_sort (const vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors, const uint64_t min_value){    
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
void newnew_bounded_counting_sort (const vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors, const uint64_t min_value){    
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


void counting_sort (const vector<uint64_t>& results, vector<pair<uint16_t, uint16_t>>& ans, const size_t found_fmin, const uint16_t n_colors){    
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

void pseudoalignment_stats(vector<int64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
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
    vector<pair<uint64_t, uint64_t>> fmin_v(fmin_counts.begin(), fmin_counts.end());
    std::sort(fmin_v.begin(), fmin_v.end()); 
    read_colors_old(data, n_colors, results, fmin_v);

    const size_t found_fmin = Fmin.size(); // # total finimizers
    counting_sort(results, ans, found_fmin, n_colors);
    return;
}

uint16_t pseudoalignment_stats(vector<int64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float& t){ 
    vector<uint64_t> tot_res;
    tot_res.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();
    
    // Count freq of each fmin
    std::unordered_map<uint64_t, uint64_t> fmin_counts;
    for (auto v : Fmin) {
        fmin_counts[v]++;
    }

    // vector for sorted output so that it is possible to scan color_set_concat

    vector<pair<uint64_t, uint64_t>> fmin_v(fmin_counts.begin(), fmin_counts.end());
    std::sort(fmin_v.begin(), fmin_v.end());

    //read_colors(data, n_colors, tot_res, fmin_v);

    // Check the values above the minimum
    const size_t found_fmin = Fmin.size(); // # total finimizers

    counting_sort(tot_res, ans, found_fmin, n_colors);
    
    // TODO use min value here?
    const uint64_t min_value = found_fmin * t;

    return min_value;
}

void foreach_union_set_bit(const uint64_t* data, uint64_t start1, uint64_t start2, uint64_t n_colors, vector<uint64_t>& results, uint64_t freq) {
    const uint64_t bit_offset1 = start1 * n_colors;
    const uint64_t bit_offset2 = start2 * n_colors;

    uint64_t idx = 0;

    while (idx < n_colors) {
        
        uint64_t word1_idx = (bit_offset1 + idx) / 64;
        uint64_t word2_idx = (bit_offset2 + idx) / 64;

        uint64_t shift1 = (bit_offset1 + idx) % 64;
        uint64_t shift2 = (bit_offset2 + idx) % 64;

        uint64_t bits_to_process = std::min({n_colors - idx, 64 - shift1, 64 - shift2});

        uint64_t mask = (bits_to_process == 64) ? ~0ULL : ((1ULL << bits_to_process) - 1);

        uint64_t word1 = (data[word1_idx] >> shift1) & mask;
        uint64_t word2 = (data[word2_idx] >> shift2) & mask;
        uint64_t combined = word1 | word2;

        while (combined) {
            uint64_t bit = __builtin_ctzll(combined);
            uint64_t result_idx = idx + bit;
            //if (result_idx < n_colors) {
                results[result_idx] += freq;
            //} 
            combined &= combined - 1;
        }

        idx += bits_to_process;
    }
}

void read_f_rc_colors(const uint64_t* data, const uint64_t n_colors, vector<uint64_t>& results, const vector<pair<pair<uint64_t, uint64_t>, uint64_t>>& p_fmin_v){
    // option 1: compare f and r as you go
    // option 2: store f and r in 2 bitvectors and compare at the end
    //for (const auto& [[start, r_start], freq] : p_fmin_v) {
    for (const auto& [key,freq] : p_fmin_v) {
        const uint64_t start = key.first;
        const uint64_t r_start = key.second;

        foreach_union_set_bit(data, start, r_start, n_colors, results, freq);
    }
}

// TODO: how to sort these or exploit identical ones still keeping the pairs?
// two overlpaiing k-mers are likely to have the same fmin so they are likely to share the same fmin on both strands
// This should anyways keep the number of false pos low
void combine_f_rc(vector<int64_t>& Fmin, vector<int64_t>& r_Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans){
    // NEW pseudoaligment_stats
    vector<uint64_t> results;
    results.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();
    // for now assume they have the same size

    //option1 
    // check the ones that are the same and deal with only the others later
    // store the diff ones in a separate vector??
    // you can still sort Fmin -> correct results
    // How to deal with the reverse??? sort r_Fmin based on Fmin sorting

    
    const int n_fmin = Fmin.size(); // now this is input_len -k +1
    // 1. Make a vector of pairs and a vector of int64_t
    vector<pair<int64_t, int64_t>> p_Fmin; // pair of distinct color_sets ids
    p_Fmin.reserve(n_fmin);

    vector<int64_t> i_Fmin; // if colorset id is the same (identical) for forward and reverse
    i_Fmin.reserve(n_fmin);

    for (auto i = 0; i<n_fmin; i++){
        int64_t f = Fmin[i];
        int64_t r = r_Fmin[n_fmin - i - 1];
        // r string is the reverse complement
        if (f != r && f != -1 && r != -1){
            p_Fmin.push_back({f, r});
        }
        else if (f != -1){
            i_Fmin.push_back(f);
        } 
        else if (r != -1){
            i_Fmin.push_back(r);
        }
    }
    // cases:
    // == -1 -1 : nothing
    // == x x : store x
    // != -1 x : store x
    // != x -1 : store x
    // != x y : store x,y

    const size_t found_fmin = p_Fmin.size() + i_Fmin.size();

    // 2. Deal with the vector of int64_t: i_Fmin
        // 2a. Keep frequency 
        // Count freq of each i_fmin
    std::unordered_map<int64_t, uint64_t> i_fmin_counts;
    i_fmin_counts.reserve(i_Fmin.size());
    for (const auto& v: i_Fmin) {
        i_fmin_counts[v]++;
    }

    // 2b. sort i_fmin
    vector<pair<int64_t, uint64_t>> i_fmin_v(i_fmin_counts.begin(), i_fmin_counts.end());
    std::sort(i_fmin_v.begin(), i_fmin_v.end()); 
    read_colors(data, n_colors, results, i_fmin_v);

    // 3. Deal with the vector of pairs: p_Fmin
    struct pair_hash {
        size_t operator()(const pair<int64_t, int64_t>& p) const {
            return std::hash<int64_t>()(p.first) ^ (std::hash<int64_t>()(p.second) << 1);
        }
    };
    std::unordered_map<pair<int64_t,int64_t>, uint64_t, pair_hash> p_fmin_counts;
    p_fmin_counts.reserve(p_Fmin.size());

    for (const auto& v: p_Fmin) {
        p_fmin_counts[v]++;
    }
    vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v(p_fmin_counts.begin(), p_fmin_counts.end()); // [{{f,r},counts},...]

    // TODO does it make sense to sort pairs??

    read_f_rc_colors(data, n_colors, results, p_fmin_v);

    counting_sort(results, ans, found_fmin, n_colors);
    return;
}

uint64_t combine_f_rc(vector<int64_t>& Fmin, vector<int64_t>& r_Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<pair<uint16_t, uint16_t>>& ans, const float& t){

    // NEW pseudoaligment_stats
    vector<uint64_t> results;
    results.resize(n_colors, 0);

    const uint64_t* data = color_sets_concat.data();


    const int n_fmin = Fmin.size(); // now this is input_len -k +1
    // 1. Make a vector of pairs and a vector of int64_t
    vector<pair<int64_t, int64_t>> p_Fmin;
    p_Fmin.reserve(n_fmin);

    vector<int64_t> i_Fmin; // if colorset id is the same for forward and reverse
    i_Fmin.reserve(n_fmin);

    for (auto i = 0; i<n_fmin; i++){
        int64_t f = Fmin[i];
        int64_t r = r_Fmin[n_fmin - i - 1];
        if (f != r && f != -1 && r != -1){
            p_Fmin.push_back({f, r});
        }
        else if (f != -1){
            i_Fmin.push_back(f);
        } 
        else if (r != -1){
            i_Fmin.push_back(r);
        }
    }
    const size_t found_fmin = p_Fmin.size() + i_Fmin.size();

    // 2. Deal with the vector of int64_t: i_Fmin
        // 2a. Keep frequency 
        // Count freq of each i_fmin
    std::unordered_map<int64_t, uint64_t> i_fmin_counts;
    i_fmin_counts.reserve(i_Fmin.size());
    for (const auto& v: i_Fmin) {
        i_fmin_counts[v]++;
    }
        // 2b. sort i_fmin
    vector<pair<int64_t, uint64_t>> i_fmin_v(i_fmin_counts.begin(), i_fmin_counts.end());
    std::sort(i_fmin_v.begin(), i_fmin_v.end()); 
    read_colors(data, n_colors, results, i_fmin_v);


    // 3. Deal with the vector of pairs: p_Fmin
    struct pair_hash {
        size_t operator()(const pair<int64_t, int64_t>& p) const {
            return std::hash<int64_t>()(p.first) ^ (std::hash<int64_t>()(p.second) << 1);
        }
    };
    std::unordered_map<pair<int64_t,int64_t>, uint64_t, pair_hash> p_fmin_counts;
    p_fmin_counts.reserve(p_Fmin.size());

    for (const auto& v: p_Fmin) {
        p_fmin_counts[v]++;
    }
    vector<pair<pair<uint64_t, uint64_t>, uint64_t>> p_fmin_v(p_fmin_counts.begin(), p_fmin_counts.end());

    // TODO does it make sense to sort pairs??
    read_f_rc_colors(data, n_colors, results, p_fmin_v);

    counting_sort(results, ans, found_fmin, n_colors);

    const uint64_t min_value = found_fmin * t;

    return min_value;
}



