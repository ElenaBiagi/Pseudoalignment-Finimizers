#pragma once

#include <string>
#include <cstring>
#include <unordered_map>
#include <limits>

#include <filesystem>
#include <cstdio>
#include <optional>
#include <deque>

//#include "PackedStrings.hh"
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

    uint64_t start = 0;
    uint64_t kmer_start = 0;
    bool found;
    uint64_t s_int;

    uint8_t f_len;
    uint64_t f_int;
    uint32_t f_color; 

    tuple<uint8_t, uint64_t, uint32_t, uint64_t> curr_substr;
    tuple<uint8_t, uint64_t, uint32_t, uint64_t> k_fmin = {k+1, 0, 0, str_len};
    BoundedDeque<tuple<uint8_t, uint64_t, uint32_t, uint64_t>> all_fmin(str_len);;
    //fmin = {f_len, f_int, f_color, start}

    /* for(auto x: sB){
        cerr << x.first << ": "<< (int)x.second.first <<endl;
    }
    cerr<< endl; */

    // iterate over input
    //cerr << str_len-plen+1<< endl;
    for (start = 0; start < str_len-plen+1; start++ ){ //TODO the last k-plen characters cannot contain a prefix
        //cerr << "START: "<< start<< endl;
        found = false;
        uint64_t int_p = prefix2int(input, start, plen);
        // 1. prefix found
        if (buckets[int_p].has_value()){
            /* cerr << "prefix found !"<< endl;
            cerr << int_p << endl;
            cerr << input.substr(start, plen) << endl;
             */
            // extract the tail
            const Bucket& bucket_p = *buckets[int_p];
            // extract the LONGEST possible tail starting from start+plen. it will be shortened by bitMagicSearch depending on tlen

            char s_len = (str_len >= start+k) ? k-plen : str_len-start-plen;
            uint64_t int_s = prefix2int(input, start+plen, s_len); // tail
            auto [pos,len] = bitMagicSearch(bucket_p.tail_data, int_s, s_len); // input: sdsl::bit_vector &T, int64_t pointer, string S    
            if (pos > -1){
                //cerr << "found! len: " << (int)len << endl;
                found = true;
                f_len = plen + len;
                f_int = (int_p << (len *2)) | (int_s >> ((s_len - len)*2) ); // TODO ADD PREFIX AND TLEN: shift p_int to the left by len*2, and int_s to the right to remove the unused chars
                f_color = bucket_p.color_set_ids[pos];
                //cerr << f_color << endl;
            }
        }else{
            //cerr << endl;
            //cerr << "prefix not found" << endl;
            // 2. Prefix NOT found
            // Start from the longest possible prefix
            // if you find a real match, stop
            uint8_t sp_len = plen-1;
            // shift int_p to the right by (plen-plen-1)*2 characters;
            //cerr << "prefix = "<< input.substr(start, plen) << " " << int_p << endl;
            uint64_t int_sp = int_p >> ((plen-sp_len)*2);
            auto it = sB.find(int_sp);
            //cerr << input.substr(start, plen-1) << endl;
            while (sp_len > 0){ // sp_len must be < plen as the whole prefix was not found
                    //cerr << "int_sp: " << int_sp << endl;

                    //cerr << "sp_len: "<< (int)sp_len << endl;
                if (it != sB.end() ){
                    
                    //cerr << "found sB" << endl;
                    //cerr << "sB len: " << (int)it->second.first << endl;
                    if (sp_len == it->second.first){ // real match
                    found = true;
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    break;
                }
            }
                sp_len--;
                int_sp = int_p >> ((plen-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }
            //cerr << endl;

        }

        if (found){
            //cerr << "found" << endl;
            //cerr << "START: "<< start<< endl;
            curr_substr = {f_len, f_int, f_color, start};
            // check later if this is a good finimizer
            all_fmin.push_back(curr_substr);
        }
        else{
            //cerr << "not found"<< endl;
        }
        if (start >= k-1){ // we have looked at all the characters of the kmer
            
            //1. check end and start are ok
            
            tuple<uint8_t, uint64_t, uint32_t, uint64_t> curr_fmin = {k+1,0,0,kmer_start};
            uint64_t i = 0;
            // if all_fmin is not empty and the while loop is not entered, k_fmin is the smallest fmin and is already ok
            // Ends after the kmer? (get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k)
            // Starts before the k-mer? get<3>(k_fmin) < kmer_start
            while ( i < all_fmin.size()){
            //loop trough all_fmin
                k_fmin = all_fmin[i]; // if this is ok it will be the smallest
                // start is now the end of the k-mer
                // len+ start of the finimizer -1 must be smaller or equal to the end of the current k-mer
                // start of finimizer must be bigger or equal start of the current k-mer 
                if ((get<0>(k_fmin) + get<3>(k_fmin) - 1 <= start) && (get<3>(k_fmin) >= kmer_start)){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                    if (curr_fmin > k_fmin){
                        curr_fmin = k_fmin;
                    }
                }
                
                i++;
            } //while loop ends here

            if (get<0>(curr_fmin) + get<3>(curr_fmin) - 1 <= start and get<3>(curr_fmin) >= kmer_start){
                Fmin.push_back(get<2>(curr_fmin));
            } // Store only that start of the colors set ids in color_set_concat
            else{
                //cerr << "kmer = " << kmer_start << ":" << start << endl;
                //cerr << "curr_fmin = " << (int)get<0>(curr_fmin)<< ", " << get<1>(curr_fmin) << ", "<< get<2>(curr_fmin) << ", "<< get<3>(curr_fmin) << endl;
                /* for (int j=0; j<all_fmin.size(); j++){
                    auto ff= all_fmin[j];
                    cerr << (int)get<0>(ff)<< ", " << get<1>(ff) << ", "<< get<2>(ff) << ", "<< get<3>(ff) << endl;
                } */
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, std::min((uint64_t)k, str_len - kmer_start)) << endl;
            }
            
            kmer_start++;
        }

    }
    //TODO check the last plen-1 values
    //cerr << "last kmers "<< endl;
    uint64_t int_p = prefix2int(input, start, plen-1);
        for (uint64_t ss = start; ss < str_len; ss++ ){ //TODO the last k-plen characters cannot contain a prefix
            //cerr << "START: "<< ss<< endl;

            // TODO add a function here
            uint8_t sp_len = 1;
            // shift int_p to the right by plen-1 characters;
            uint64_t int_sp = int_p >> ((plen-1-sp_len)*2);
            auto it = sB.find(int_sp);

            while (sp_len < plen-1){ // sp_len must be < plen as the whole prefix was not found
                if (it != sB.end() && sp_len == it->second.first){ // real match
                    found = true;
                    f_len = sp_len;
                    f_int = int_sp;
                    f_color = it->second.second;
                    break;
                }
                sp_len++;
                int_sp = int_p >> ((plen-1-sp_len)*2); // shift to the right by 2 character less than then prev time
                it = sB.find(int_sp);
            }

            if (found){
                //cerr << "found!"<< endl;
                curr_substr = {f_len, f_int, f_color, ss};
                all_fmin.push_back(curr_substr);
            }
            if (ss >= k-1){ // we have looked at all the characters of the kmer
            //1. check end is ok
            // TODO is it not ok to check end only before? no becuase we are modifying k_fmin here
            tuple<uint8_t, uint64_t, uint32_t, uint64_t> curr_fmin = {k+1,0,0,kmer_start};
            uint64_t i = 0;
            // if all_fmin is not empty and the while loop is not entered, k_fmin is the smallest fmin and is already ok
            // Ends after the kmer? (get<0>(k_fmin) + get<3>(k_fmin) >= kmer_start+k)
            // Starts before the k-mer? get<3>(k_fmin) < kmer_start
            while ( i < all_fmin.size()){
            //loop trough all_fmin
                k_fmin = all_fmin[i]; // if this is ok it will be the smallest
                // start is now the end of the k-mer
                // len+ start of the finimizer -1 must be smaller or equal to the end of the current k-mer
                // start of finimizer must be bigger or equal start of the current k-mer 
                if ((get<0>(k_fmin) + get<3>(k_fmin) - 1 <= ss) && get<3>(k_fmin) >= kmer_start){ // {length, fmin, C_offset, start} // if start comes before the kmer_start that it must be discarded
                    if (curr_fmin > k_fmin){
                        curr_fmin = k_fmin;
                    }
                }
                
                i++;
            }
            if (get<0>(curr_fmin) + get<3>(curr_fmin) - 1 <= ss and get<3>(curr_fmin) >= kmer_start){
                Fmin.push_back(get<2>(curr_fmin));
            }
            else{
                /* for (int j=0; j<all_fmin.size(); j++){
                    auto ff= all_fmin[j];
                    //cerr << (int)get<0>(ff)<< ", " << get<1>(ff) << ", "<< get<2>(ff) << ", "<< get<3>(ff) << endl;
                } */
                cerr << "finimizer not found for kmer " << kmer_start << " " << input.substr(kmer_start, std::min((uint64_t)k, str_len - kmer_start)) << endl;}
            
            //2. check that start is ok for the next kmer
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
        cerr << Fmin.size() << endl;
        results.assign(n_colors, 0);
        
        //#pragma omp parallel for
        for(const auto& start : Fmin){
            cerr << start << endl;
            for (uint64_t i=0; i< n_colors; i++){
                //#pragma omp atomic
                results[i]+=color_sets_concat[start+i];
                //cerr << color_sets_concat[start+i];
            }
        }
        //cerr << endl;
        return;
    }

    void pseudoalignemnt_stats(const vector<uint64_t>& Fmin, const sdsl::bit_vector& color_sets_concat, const uint64_t n_colors, vector<float>& results, const float& t){ 
        // count the number of finimizers found
        size_t found_fmin = Fmin.size(); // # total finimizers
        //size_t rm_fmin = 0; // finimizers not found in the index
        //cerr << Fmin.size() << endl;
        
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

//TODO rewrite
// used in build
/* string print_finimizer_stats(const set<tuple<int64_t, int64_t, int64_t>>& finimizers, int64_t n_kmers, int64_t n_nodes, int64_t t){
    int64_t new_number_of_fmin = finimizers.size();
    int64_t sum_freq = 0;
    int64_t sum_len = 0;
    for (auto x : finimizers){
        sum_freq += get<1>(x);
        sum_len += get<0>(x);
    }

    string result = to_string(new_number_of_fmin) + "," + to_string(sum_freq) + "," + to_string(static_cast<float>(sum_freq) / static_cast<float>(new_number_of_fmin)) + "," + to_string(static_cast<float>(sum_len) / static_cast<float>(new_number_of_fmin)) + "," + to_string(n_kmers);

    write_log(to_string(t) + "," + result, LogLevel::MAJOR);
    write_log("#Distinct finimizers: " + to_string(new_number_of_fmin) , LogLevel::MAJOR);
    write_log("Sum of frequencies: " + to_string(sum_freq) , LogLevel::MAJOR);
    write_log("Avg frequency: " + to_string(static_cast<float>(sum_freq)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    write_log("Avg length: " + to_string(static_cast<float>(sum_len)/static_cast<float>(new_number_of_fmin)) , LogLevel::MAJOR);
    return result;
} */