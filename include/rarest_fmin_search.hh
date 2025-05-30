#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

#include <omp.h>


#include "SeqIO.hh"
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

void FindPrefix(const vector<optional<Bucket>>& buckets, const uint8_t plen, const char s_len, const uint64_t int_s, const uint64_t int_p, const uint64_t start, CBuffer& all_fmin ){
    const Bucket& bucket_p = *buckets[int_p];
    auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
    if (pos > -1){
        uint64_t f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
        all_fmin.insert(make_tuple(plen + len, f_int, bucket_p.color_set_ids[pos], start));
    }
}
            
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 
    const int64_t str_len = input.size();
    if (str_len < k){return {};}
    vector<uint64_t> Fmin;// pointer to C

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

    if (buckets[int_p].has_value()){
        // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen
        // if the length of input is at leas tk , the slen = k-plen
        // defined above int_s = prefix2int(input, start+plen, s_len); // tail
        FindPrefix(buckets, plen, s_len, int_s, int_p, start, all_fmin); 
    }else{
        int_sp = int_p >> 2;
        FindShortFinimizer(plen-1, int_sp, sB, start, all_fmin);
    }

    // iterate over input
    for (start = 1; start < str_len-plen+1; start++ ){ //TODO the last k-plen characters cannot contain a prefix
        // 1. prefix found
        int_p = roll_kmer(int_p, input[start + plen - 1], plen);
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
  
        if (start >= k-1){ // we have looked at all the characters of the kmer
            
            //1. check end and start are ok
            MyTuple best_fmin = make_tuple(k+1,0,0,kmer_start);
            
            all_fmin.for_each_recent([&start, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
                // start of finimizer must be bigger or equal start of the current k-mer 
                const auto& [f_len, f_int, f_color, f_start] = k_fmin;
                if ((f_start >= kmer_start) && (f_len + f_start - 1 <= start) && (best_fmin > k_fmin) ){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                        best_fmin = k_fmin;
                    } 
                    /* return true;
                } else{ return false; }  */
            });

            if (get<0>(best_fmin) < k+1){Fmin.push_back(get<2>(best_fmin));} // Store only the start of the color set ids in color_set_concat
            else{
                /* all_fmin.for_each_recent([&start, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
                // start of finimizer must be bigger or equal start of the current k-mer 
                const auto& [f_len, f_int, f_color, f_start] = k_fmin;
                cerr << "{"<< (int)f_len << ", "<< f_start << " }" << endl;
            }); */
                // TODO remove
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, min((uint64_t)k, str_len - kmer_start)) << endl;
            }
            
            kmer_start++;
        }

    }
    // Check the last plen-1 values
    uint8_t s_plen = plen;
    for (uint64_t ss = start; ss < str_len; ss++ ){ //TODO the last k-plen characters cannot contain a prefix
        s_plen--;
        int_p = prefix2int(input, ss, s_plen); // shorten by 1 at every loop iteration
        FindShortFinimizer(s_plen, int_p, sB, ss, all_fmin);// this shortens s_plen by 1 internally
        
        if (ss >= k-1){ // we have looked at all the characters of the kmer
            //1. check end is ok
            MyTuple best_fmin = make_tuple(k+1,0,0,kmer_start);
            all_fmin.for_each_recent([&ss, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
                const auto& [f_len, f_int, f_color, f_start] = k_fmin;
                if ((f_start >= kmer_start) && (f_len + f_start - 1 <= ss) && (best_fmin > k_fmin) ){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                        best_fmin = k_fmin;
                    } 
                /*     return true;
                } else{ return false; } */ 
            });
            if (get<0>(best_fmin) < k+1){Fmin.push_back(get<2>(best_fmin));} // Store only the start of the color set ids in color_set_concat
            else{
                // TODO remove
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, min((uint64_t)k, str_len - kmer_start)) << endl;
            }
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
        
        vector<vector<uint64_t>> local_results(omp_get_max_threads(), vector<uint64_t>(n_colors, 0));

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto& local = local_results[tid];

            #pragma omp for
            for (size_t j = 0; j < Fmin.size(); ++j) {
                uint64_t start = Fmin[j];
                uint64_t base = n_colors * start;
                for (uint64_t i = 0; i < n_colors; ++i) {
                    local[i] += color_sets_concat[base + i];
                }
            }
        }

        // sum final results
        for (int t = 0; t < local_results.size(); ++t) {
            for (uint64_t i = 0; i < n_colors; ++i) {
                results[i] += local_results[t][i];
            }
        }
    }


    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
        size_t found_fmin = Fmin.size();
        vector<vector<uint64_t>> local_results(omp_get_max_threads(), vector<uint64_t>(n_colors, 0));

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto& local = local_results[tid];

            #pragma omp for
            for (size_t j = 0; j < Fmin.size(); ++j) {
                uint64_t start = Fmin[j];
                for (uint64_t i = 0; i < n_colors; ++i) {
                    local[i] += color_sets_concat[start + i];
                }
            }
        }

        vector<uint64_t> tot_res(n_colors, 0);
        for (const auto& local : local_results) {
            for (uint64_t i = 0; i < n_colors; ++i) {
                tot_res[i] += local[i];
            }
        }

        results.resize(n_colors);
        #pragma omp parallel for
        for (uint64_t i = 0; i < n_colors; ++i) {
            float fraction = static_cast<float>(tot_res[i]) / found_fmin;
            results[i] = (t == 1.0f && fraction == 1.0f) ? 1.0f : ((fraction > t) ? fraction : 0.0f);
        }
    }