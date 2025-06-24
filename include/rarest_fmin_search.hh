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
            if ((start + sp_len -1) > end) { next_candidates.push_back(make_tuple(sp_len+start-1, int_sp, it->second.second, start)); } // Sorted based on END
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
    auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
    if (pos > -1){
        uint64_t f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        //all_fmin.insert(make_tuple(plen + len, f_int, bucket_p.color_set_ids[pos], start));
        if ((start + plen+len -1) > end) { next_candidates.push_back(make_tuple(plen+len+start-1, f_int, bucket_p.color_set_ids[pos], start)); } // Sorted based on END
        else { 
            tuple<uint64_t, uint64_t, uint64_t, uint64_t> new_fmin = {plen + len, f_int, bucket_p.color_set_ids[pos], start};
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
    }
}

void PickFinimizer(vector<uint64_t>& Fmin, const uint64_t kmer_start, const uint64_t k, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& curr_candidates, BoundedDeque<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>& next_candidates, tuple<uint64_t, uint64_t, uint64_t, uint64_t>& k_fmin){
    if (!curr_candidates.empty()){
        k_fmin= curr_candidates.front();
        Fmin.push_back(get<2>(k_fmin));
    }
/*     else{
        // TODO remove this branch
        cerr << "finimizer not found for kmer " << kmer_start << endl;//<< " " << input.substr(kmer_start, min((uint64_t)k, str_len - kmer_start)) << endl;
    } */

    // 1. Check if this finimizer is good for the next k-mer (still in the window)
    if (get<3>(k_fmin) == kmer_start){ // this is never true if curr_size is empty
        curr_candidates.pop_front();
        k_fmin = (curr_candidates.empty()) ? static_cast<tuple<uint64_t, uint64_t, uint64_t, uint64_t>>(make_tuple(k+1,0,0,kmer_start+1)) : curr_candidates.front();
    }

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
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint64_t plen, const uint64_t k){ 

    const int64_t str_len = input.size();
    //if (str_len < k){return {};} // This is checked before already
    vector<uint64_t> Fmin;// pointer to C

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
        FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); 
    }else{
        // 2. look for a shorter finimizer
        int_sp = int_p >> 2;
        FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); 
    }
    
    // The first k-1 characters do not contail all possible finimizers for the first k-mer
    for (start = 1; start < k-1; start++ ){ 
        int_p = stream_kmer(int_p, input[start + plen - 1], plen); // shorten by 1 at every loop iteration 
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
        int_s = prefix2int(input, start+plen, s_len);; // tail   
        if (buckets[int_p].has_value()){
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); 
        }
        else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); 
        }
    }

    uint64_t ss = start;
    s_len = k - plen; // Constant s_len
    for (start = ss ; start <= str_len-k; start++ ){ 
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen 
        int_s = stream_kmer(int_s, input[start + plen + s_len - 1], s_len);
        
        if (buckets[int_p].has_value()){
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin);  
        }else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin);  
        }
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); 
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
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin);  
        }else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin);  
        }
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); 
        kmer_start++;
    }

    //The last plen-1 characters cannot contain a prefix
    ss = start;
    uint64_t s_plen = plen;
    for (start = ss; start < str_len; start++ ){ 
        s_plen--;
        int_p &= ((1ULL << (2 * s_plen)) - 1); // Shorten int_p by 2 
        FindShortFinimizer(s_plen, int_p, sB, start, kmer_start+k-1, curr_candidates, next_candidates, k_fmin); // this shortens s_plen by 1 internally
        
        PickFinimizer(Fmin, kmer_start, k, curr_candidates, next_candidates, k_fmin); 
        kmer_start++;
    }
    
    return Fmin;
}


void pseudoalignment_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<uint64_t>& results){ 
    results.assign(n_colors, 0);
    const uint64_t* data = color_sets_concat.data();

    for(const auto& start : Fmin){
        size_t offset_bits = n_colors*start;
        size_t word_index = offset_bits / 64;
        size_t w_offset = offset_bits % 64;
        const uint64_t bits_read = std::min(64 - w_offset, n_colors);

        uint64_t color_id = 0;

        // 1. Read the first word
        uint64_t mask = (bits_read == 64) ? ~0ULL : ((1ULL << bits_read) - 1);
        uint64_t first_w = (data[word_index] >> w_offset) & mask;
        
        while (first_w != 0) {
            uint64_t lowest_set_bit = __builtin_ctzll(first_w);
            results[lowest_set_bit]++;
            first_w &= first_w - 1;
        }
        color_id += bits_read;


        // 2. Read aligned words in btw
        w_offset = 0;
        // Check how many words we will have to read
        uint64_t bits_left = n_colors - color_id;
        while (color_id + 64 <= n_colors) {
            word_index++;
            uint64_t w = data[word_index];
            while (w != 0) {
                uint64_t lowest_set_bit = __builtin_ctzll(w);
                results[color_id + lowest_set_bit]++;
                w &= w - 1;
            }
            color_id += 64;
        }

        // 3. Read the last word (if any)
        bits_left = n_colors - color_id;
        if (bits_left > 0){
            uint64_t mask = (bits_left == 64) ? ~0ULL : ((1ULL << bits_left) - 1);
            uint64_t last_w = data[word_index + 1] & mask;

            while (last_w != 0) {
                uint64_t lowest_set_bit = __builtin_ctzll(last_w);
                results[color_id + lowest_set_bit]++;
                last_w &= last_w - 1;
            }
        }
    }
    return;
}

void pseudoalignment_stats(vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
    std::sort(Fmin.begin(), Fmin.end());    
    size_t found_fmin = Fmin.size(); // # total finimizers

    vector<uint64_t> tot_res;
    tot_res.assign(n_colors, 0);

    for(const auto& start : Fmin){
        for (uint64_t i=0; i< n_colors; i++){
            tot_res[i]+=color_sets_concat[start+i];
        }
    }
    results.assign(n_colors, 0);

    for (uint64_t i=0; i< n_colors; i++){         
        // For every color found, (#finimizers with that color)/(#tot finimizers)
        float fraction = static_cast<float>(tot_res[i]/static_cast<float>(found_fmin));

            if (fraction >= t){
                results[i]=fraction;
            }
        }
        return;
    }
 