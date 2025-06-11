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
//#include "BoundedDeque.hh"
#include "CircularBuffer.hh"

#include "common.hh"
#include "bitsearch.hh"
#include "Buckets.hh"

using MyTuple = tuple<uint8_t, uint64_t, uint32_t, uint64_t>; // {f_len, f_int, f_color, start}


void FindShortFinimizer(const uint8_t int_sp_len, uint64_t int_sp, const unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint64_t start, CBuffer& all_fmin){
    // 2. Prefix NOT found
    // Start from the longest possible prefix
    // if you find a real match, stop
    uint8_t sp_len = int_sp_len; // plen is here plen-1: sp_len must be < plen as the whole prefix was not found  
    while (sp_len > 0){ 
        auto it = sB.find(int_sp);
        if (it != sB.end() && sp_len == it->second.first){ // real match 
            all_fmin.insert(make_tuple(sp_len, int_sp, it->second.second, start));
            return;
        }
        int_sp >>= 2;
        sp_len--;
    }
}

void FindPrefix(const vector<optional<Bucket>>& buckets, const uint8_t plen, const char s_len, const uint64_t int_s, const uint64_t int_p, const uint64_t start, CBuffer& all_fmin){
    const Bucket& bucket_p = *buckets[int_p];
    auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
    if (pos > -1){
        uint64_t f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        all_fmin.insert(make_tuple(plen + len, f_int, bucket_p.color_set_ids[pos], start));
    }
}

void PickFinimizer(vector<uint64_t>& Fmin, const uint64_t kmer_start, CBuffer& all_fmin,const uint8_t k, vector<tuple<int64_t, int64_t, int64_t>>& finimizers){
    MyTuple best_fmin = make_tuple(k+1,0,0,kmer_start);
    all_fmin.for_each_recent([&k, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
    // start of finimizer must be bigger or equal start of the current k-mer 
    const auto& [f_len, f_int, f_color, f_start] = k_fmin;
    if ((f_start >= kmer_start) && (f_len + f_start - 1 <= (kmer_start+k-1)) && (best_fmin > k_fmin) ){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
            best_fmin = k_fmin;
        } 
        /* return true;
    } else{ return false; }  */
    });

    if (get<0>(best_fmin) < k+1){
        Fmin.push_back(get<2>(best_fmin));

        tuple<int64_t, int64_t, int64_t> curr_fmin = {get<0>(best_fmin),get<1>(best_fmin),get<2>(best_fmin)};

        if (!finimizers.empty()){
            if (curr_fmin != finimizers.back()){
                finimizers.push_back(curr_fmin);
            } // Store only the start of the color set ids in color_set_concat
        } 
        else {
            finimizers.push_back(curr_fmin);
        }
    }
    else{
        // TODO remove this branch
        /* all_fmin.for_each_recent([&start, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
        // start of finimizer must be bigger or equal start of the current k-mer 
        const auto& [f_len, f_int, f_color, f_start] = k_fmin;
        cerr << "{"<< (int)f_len << ", "<< f_start << " }" << endl;
        }); */
        cerr << "finimizer not found for kmer "<< endl;// << kmer_start << " " << input.substr(kmer_start, min((uint64_t)k, str_len - kmer_start)) << endl;
    }
}
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 

    const int64_t str_len = input.size();
    if (str_len < k){return {};}
    vector<uint64_t> Fmin;// pointer to C
    vector<tuple<int64_t, int64_t, int64_t>> finimizers;

    uint64_t start = 0;
    uint64_t kmer_start = 0;
    uint64_t s_int;

    uint8_t f_len;
    uint64_t f_int;
    uint32_t f_color; 

    MyTuple curr_substr;
    MyTuple k_fmin = make_tuple(k+1, 0, 0, str_len);
    //vector<MyTuple> all_fmin(str_len);
    CBuffer all_fmin(k);
    
    uint64_t int_p = prefix2int(input, start, plen); // start = 0
    char s_len = k-plen;
    uint64_t int_s = prefix2int(input, start+plen, s_len); // tail
    uint64_t int_sp;

    // Check the first k-2 characters
    if (buckets[int_p].has_value()){
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        // if the length of input is at leas k, slen = k-plen
        // defined above int_s = prefix2int(input, start+plen, s_len); // tail
        FindPrefix(buckets, plen, s_len, int_s, int_p, start, all_fmin); 
    }else{
        int_sp = int_p >> 2;
        FindShortFinimizer(plen-1, int_sp, sB, start, all_fmin);
    }
    
     
    // The first k-1 characters do not contail all possible finimizers for the first k-mer
    for (start = 1; start < k-1; start++ ){ 
        // 1. prefix found
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        //int_p = prefix2int(input, start, plen); // shorten by 1 at every loop iteration    
        if (buckets[int_p].has_value()){
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
            s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
            int_s = prefix2int(input, start+plen, s_len);; // tail
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, all_fmin);
        }else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, all_fmin);
        }
    }
    uint64_t ss = start;
    for (start = ss ; start < str_len-plen+1; start++ ){ // TODO end this loop earlier so that s_len is constant 
        // 1. prefix found
        int_p = stream_kmer(int_p, input[start + plen - 1], plen);
        //int_p = prefix2int(input, start, plen); // shorten by 1 at every loop iteration    
        if (buckets[int_p].has_value()){
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
            s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
            int_s = prefix2int(input, start+plen, s_len);; // tail
            FindPrefix(buckets, plen, s_len, int_s, int_p, start, all_fmin);
        }else{
            int_sp = int_p >> 2;
            FindShortFinimizer(plen-1, int_sp, sB, start, all_fmin);
        }
        PickFinimizer(Fmin, kmer_start, all_fmin, k, finimizers);
        kmer_start++;
    }
    // TODO Check the lenght of the string and add another loop for shorter s_len and stream s_int

    // Check the last plen-1 values
    //The last k-plen characters cannot contain a prefix
    ss = start;
    uint8_t s_plen = plen;
    for (start = ss; start < str_len; start++ ){ 
        s_plen--;
        int_p = prefix2int(input, start, s_plen); // shorten by 1 at every loop iteration
        FindShortFinimizer(s_plen, int_p, sB, start, all_fmin);// this shortens s_plen by 1 internally
        
        if (start >= k-1){ // we have looked at all the characters of the kmer
            PickFinimizer(Fmin, kmer_start, all_fmin, k, finimizers);
            kmer_start++;
        }
    }
    
    return Fmin;
}


void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<uint64_t>& results){ 
    
    results.assign(n_colors, 0);
    const uint64_t* data = color_sets_concat.data();

    for(const auto& start : Fmin){
        size_t offset_bits = n_colors*start;
        size_t word_index = offset_bits / 64;
        size_t w_offset = offset_bits % 64;
        const uint8_t bits_read = std::min(64 - w_offset, n_colors);

        uint64_t color_id = 0;

        // 1. Read the first word
        uint64_t mask = (bits_read == 64) ? ~0ULL : ((1ULL << bits_read) - 1);
        uint64_t first_w = (data[word_index] >> w_offset) & mask;
        
        while (first_w != 0) {
            uint8_t lowest_set_bit = __builtin_ctzll(first_w);
            results[lowest_set_bit]++;
            first_w &= first_w - 1;
        }
        color_id += bits_read;


        // 2. Read aligned words in btw
        w_offset = 0;
        // Check how many words we will have to read
        uint64_t bits_left = n_colors - color_id;
        while (color_id + 64 <= n_colors - bits_left) {
            word_index++;
            uint64_t w = data[word_index];
            while (w != 0) {
                uint8_t lowest_set_bit = __builtin_ctzll(w);
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
                uint8_t lowest_set_bit = __builtin_ctzll(last_w);
                results[color_id + lowest_set_bit]++;
                last_w &= last_w - 1;
            }
        }
    }
    return;
}

void pseudoalignemnt_stats(vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
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
 