#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

#include "SeqIO.hh"
//#include "BoundedDeque.hh"
#include "CircularBuffer.hh"

#include "common.hh"
#include "bitsearch.hh"
#include "Buckets.hh"

using MyTuple = std::tuple<uint8_t, uint64_t, uint32_t, uint64_t>; // {f_len, f_int, f_color, start}


// TODO simplify this removing what is not necessary
// Do we want to count the number of found kmers? YES
// set ?
vector<uint64_t> rarest_fmin_streaming_search(const string& input, const vector<optional<Bucket>>& buckets, const std::unordered_map<uint32_t, pair<char, int64_t>>& sB, const uint8_t plen, const uint8_t k){ 
    const int64_t str_len = input.size();
    vector<uint64_t> Fmin;// pointer to C

    uint64_t start = 0;
    uint64_t kmer_start = 0;
    uint64_t s_int;

    uint8_t f_len;
    uint64_t f_int;
    uint32_t f_color; 

    MyTuple curr_substr;
    MyTuple k_fmin = std::make_tuple(k+1, 0, 0, str_len);
    //vector<MyTuple> all_fmin(str_len);
    CBuffer all_fmin(k);
    

    // iterate over input
    for (start = 0; start < str_len-plen+1; start++ ){ //TODO the last k-plen characters cannot contain a prefix
        uint64_t int_p = prefix2int(input, start, plen);
        // 1. prefix found
        if (buckets[int_p].has_value()){

            // extract the tail
            const Bucket& bucket_p = *buckets[int_p];
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen

            char s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
            uint64_t int_s = prefix2int(input, start+plen, s_len); // tail
            auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            if (pos > -1){
                f_len = plen + len;
                f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
                f_color = bucket_p.color_set_ids[pos];
                all_fmin.insert(std::make_tuple(f_len, f_int, f_color, start));

            }
        }else{
            // 2. Prefix NOT found
            // Start from the longest possible prefix
            // if you find a real match, stop
            uint8_t sp_len = plen-1;
            // shift int_p to the right by (plen-plen-1)*2 characters;
            uint64_t int_sp = int_p >> ((plen-sp_len)*2);
            auto it = sB.find(int_sp);
            while (sp_len > 0){ // sp_len must be < plen as the whole prefix was not found
                    
                if (it != sB.end() ){
                    
                    if (sp_len == it->second.first){ // real match
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    all_fmin.insert(std::make_tuple(f_len, f_int, f_color, start));
                    break;
                }
            }
                sp_len--;
                int_sp = int_p >> ((plen-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }

        }
  
        if (start >= k-1){ // we have looked at all the characters of the kmer
            
            //1. check end and start are ok
            
            MyTuple best_fmin = std::make_tuple(k+1,0,0,kmer_start);
            // if all_fmin is not empty and the while loop is not entered, k_fmin is the smallest fmin and is already ok
            // Ends after the kmer? (get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k)
            // Starts before the k-mer? get<3>(k_fmin) < kmer_start
            // TODO IMPROVE THIS
            // we could scann only the last k positions
            // we could store only the last k pos
            // LINKED LIST?? i do not want random access. I want to look at each of them
            // use a circular buffer
            
            all_fmin.for_each_recent([&start, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
                // start of finimizer must be bigger or equal start of the current k-mer 
                if ((get<3>(k_fmin) >= kmer_start)){
                    // start is now the end of the k-mer
                    // len+ start of the finimizer -1 must be smaller or equal to the end of the current k-mer
                    if ((get<0>(k_fmin) + get<3>(k_fmin) - 1 <= start) && (best_fmin > k_fmin) ){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                        best_fmin = k_fmin;
                    } 
                    return true;
                } else{ return false; } 
            });

            if (get<0>(best_fmin) < k+1){Fmin.push_back(get<2>(best_fmin));} // Store only the start of the color set ids in color_set_concat
            /* else{
                // TODO remove
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, std::min((uint64_t)k, str_len - kmer_start)) << endl;
            } */
            
            kmer_start++;
        }

    }
    // Check the last plen-1 values
    uint64_t int_p = prefix2int(input, start, plen-1);
    for (uint64_t ss = start; ss < str_len; ss++ ){ //TODO the last k-plen characters cannot contain a prefix
        
        // TODO add a function here
        uint8_t sp_len = 1;
        // shift int_p to the right by plen-1 characters;
        uint64_t int_sp = int_p >> ((plen-1-sp_len)*2);
        auto it = sB.find(int_sp);

        while (sp_len < plen-1){ // sp_len must be < plen as the whole prefix was not found
            if (it != sB.end() && sp_len == it->second.first){ // real match
                f_len = sp_len;
                f_int = int_sp;
                f_color = it->second.second;
                all_fmin.insert(std::make_tuple(f_len, f_int, f_color, ss));
                break;
            }
            sp_len++;
            int_sp = int_p >> ((plen-1-sp_len)*2); // shift to the right by 2 character less than then prev time
            it = sB.find(int_sp);
        }

        
        if (ss >= k-1){ // we have looked at all the characters of the kmer
            //1. check end is ok
            MyTuple best_fmin = std::make_tuple(k+1,0,0,kmer_start);
            // if all_fmin is not empty and the while loop is not entered, k_fmin is the smallest fmin and is already ok
            // Ends after the kmer? (get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k)
            // Starts before the k-mer? get<3>(k_fmin) < kmer_start
            // TODO IMPROVE THIS
            //for (auto k_fmin : all_fmin){
            all_fmin.for_each_recent([&ss, &kmer_start, &best_fmin](const MyTuple& k_fmin) {            
                // start of finimizer must be bigger or equal start of the current k-mer 
                if ((get<3>(k_fmin) >= kmer_start)){
                    // start is now the end of the k-mer
                    // len+ start of the finimizer -1 must be smaller or equal to the end of the current k-mer
                    if ((get<0>(k_fmin) + get<3>(k_fmin) - 1 <= ss) && (best_fmin > k_fmin) ){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                        best_fmin = k_fmin;
                    } 
                    return true;
                } else{ return false; } 
            });
            
            if (get<0>(best_fmin) < k+1){Fmin.push_back(get<2>(best_fmin));} // Store only the start of the color set ids in color_set_concat
            /* else{
                // TODO remove
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, std::min((uint64_t)k, str_len - kmer_start)) << endl;
            } */
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
                results[i]+=color_sets_concat[(n_colors*start)+i];
            }
        }
        
        return;
    }

    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
        size_t found_fmin = Fmin.size(); // # total finimizers

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